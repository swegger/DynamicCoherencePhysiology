function predictDynCohRRdynamicsModel(subject,varargin)
%%
%
%
%
%%

%% Defaults

averagingWin = [600,800];
correlationWin = [600,800];


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'neuronTypingFile',NaN)
addParameter(Parser,'fitReducedRankDynamicsFile',NaN)
addParameter(Parser,'MTtoFEFregressionFile',NaN)
addParameter(Parser,'speedInd',2)
addParameter(Parser,'cohInds',[2 1 3])
addParameter(Parser,'transitionTime',490)
addParameter(Parser,'timeconstant',20)
addParameter(Parser,'inputGain',[8 2])
addParameter(Parser,'useMTdata',true)
addParameter(Parser,'cohs',[20,60,100])
addParameter(Parser,'deltaCoh_exp',[0,-40,40])
addParameter(Parser,'saveResults',false)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'mirrorPositivePerturbation',false)
addParameter(Parser,'motionChangeTimes',[150+0:300:1500])
addParameter(Parser,'cohLookUpTable',[5,4,3])
addParameter(Parser,'dynDataLoad','neuronTyping')
addParameter(Parser,'trialCutoff',0)
addParameter(Parser,'combo',{[5],[2,4],[1,3]});
addParameter(Parser,'averagingWin',[600,800])
addParameter(Parser,'correlationWin',[600,800])

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
neuronTypingFile = Parser.Results.neuronTypingFile;
fitReducedRankDynamicsFile = Parser.Results.fitReducedRankDynamicsFile;
MTtoFEFregressionFile = Parser.Results.MTtoFEFregressionFile;
speedInd = Parser.Results.speedInd;
cohInds = Parser.Results.cohInds;
transitionTime = Parser.Results.transitionTime;
timeconstant = Parser.Results.timeconstant;
inputGain = Parser.Results.inputGain;
useMTdata = Parser.Results.useMTdata;
cohs = Parser.Results.cohs;
deltaCoh_exp= Parser.Results.deltaCoh_exp;
saveResults = Parser.Results.saveResults;
saveFigures = Parser.Results.saveFigures;
mirrorPositivePerturbation = Parser.Results.mirrorPositivePerturbation;
motionChangeTimes = Parser.Results.motionChangeTimes;
cohLookUpTable = Parser.Results.cohLookUpTable;
dynDataLoad = Parser.Results.dynDataLoad;
trialCutoff = Parser.Results.trialCutoff;
combo = Parser.Results.combo;
averagingWin = Parser.Results.averagingWin;
correlationWin = Parser.Results.correlationWin;

%% Check for specific files as input, otherwise set default files
if any(isnan(neuronTypingFile))
    switch subject
        case 'ar'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
            predictFiringRatesFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/predictFiringRatesRRdynamics20240819.mat';
        case 'fr'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';            
            predictFiringRatesFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/predictFiringRatesRRdynamics20240822.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(fitReducedRankDynamicsFile))
    switch subject
        case 'ar'
%             fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240613.mat';
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240729.mat';
        case 'fr'
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/fitReducedRankModel20240613.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(MTtoFEFregressionFile))
    switch subject
        case 'ar'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/ar/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        case 'fr'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/fr/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

%% Load the data
neuronTyping = load(neuronTypingFile,'subjects','Y','cellID','gainRegression','idx','Rinit','NumClusters','initCoh','Rdyn');
predictFiringRates = load(predictFiringRatesFile,'B','cc');
rrModel = load(fitReducedRankDynamicsFile,'modelFEF','theoreticalInput','kappas');
MTtoFEF = load(MTtoFEFregressionFile,'vOpt');
t = rrModel.modelFEF.t;
modelFEF = rrModel.modelFEF;

%% Create weighted MT response for 60% coh to transitions to 20 and 100% coh
if useMTdata
    % Model the change in inputs as the smooth change from before coherence
    % change to after, plus the extra response of MT neurons weighted by
    % log2(preferred speed) or overall sum.
    [weightedAverage, standardAverage, interpolatedWeightedAverage, interpolatedStandardAverage, conditions, mt_t, figureHandles] = analyzeMTcoherencePertubationData();
    
    dweightedAverage = weightedAverage - interpolatedWeightedAverage;
    dstandardAverage = standardAverage - interpolatedStandardAverage;
    
    pertInds = conditions(:,3) ~= 0;
    speedsMT = conditions(pertInds,1);
    deltaCohs = conditions(pertInds,3);
    
    [dCs,Ss] = meshgrid(unique(deltaCohs),unique(speedsMT));
%     Ss = Ss';
%     dCs = dCs';
    [dCsq,Ssq] = meshgrid(deltaCoh_exp(2:end),10);
%     Ssq = Ssq';
%     dCsq = dCsq';
    
    
    for ti = 1:size(weightedAverage,1)
        for i = 1:size(Ss,1)
            for j = 1:size(Ss,2)
                dWA(i,j) = dweightedAverage(ti,conditions(:,1) == Ss(i,j) & conditions(:,3) == dCs(i,j));
            end
        end
        B = regress(dWA(:),[Ss(:),dCs(:),ones(size(Ss(:)))]);
        est = [Ssq(:), dCsq(:), ones(size(Ssq(:)))]*B;
        pertR(ti,:,:,1) = reshape(est,size(Ssq));
        
        for i = 1:size(Ss,1)
            for j = 1:size(Ss,2)
                dSA(i,j) = dstandardAverage(ti,conditions(:,1) == Ss(i,j) & conditions(:,3) == dCs(i,j));
            end
        end
        B = regress(dSA(:),[Ss(:),dCs(:),ones(size(Ss(:)))]);
        est = [Ssq(:), dCsq(:), ones(size(Ssq(:)))]*B;
        pertR(ti,:,:,2) = reshape(est,size(Ssq));
        
    end
    
    % Interpolate between the two conditions
    inputInitCoh = rrModel.theoreticalInput;
    input = inputInitCoh(:,:,speedInd,cohInds(1));
    transitionInd = find(t == transitionTime);
    for i = 2:length(cohInds)
        input(1,:,i) = inputStreamInterpolate2(inputInitCoh(1,:,speedInd,cohInds(1)),...
            inputInitCoh(1,:,speedInd,cohInds(i)),...
            transitionInd,timeconstant);
        input(2,:,i) = inputStreamInterpolate2(inputInitCoh(2,:,speedInd,cohInds(1)),...
            inputInitCoh(2,:,speedInd,cohInds(i)),...
            transitionInd,timeconstant);
    end
    
    pertRTemp = pertR;
    if mirrorPositivePerturbation
        pertRTemp(:,:,1,:) = -pertR(:,:,2,:);
    end
    
    % Add the perturbation
    for i = 2:length(cohInds)
        input(1,t>=transitionTime & t<=transitionTime+(mt_t(end)-500),i) = input(1,t>=transitionTime & t<=transitionTime+(mt_t(end)-500),i) + ...
            interp1(mt_t(mt_t>=500),pertRTemp(mt_t>=500,1,i-1,1),500:mt_t(end));
        input(2,t>=transitionTime & t<=transitionTime+(mt_t(end)-500),i) = input(2,t>=transitionTime & t<=transitionTime+(mt_t(end)-500),i) + ...
            interp1(mt_t(mt_t>=500),pertRTemp(mt_t>=500,1,i-1,2),500:mt_t(end));
    end
    input(3,:,:) = inputInitCoh(3,:,3,1:3);
    
else
    inputInitCoh = rrModel.theoreticalInput;
    input = inputInitCoh(:,:,speedInd,cohInds(1));
    transitionInd = find(t == transitionTime);
    for i = 2:length(cohInds)
        input(1,:,i) = inputStreamInterpolate(inputInitCoh(1,:,speedInd,cohInds(1)),...
            inputInitCoh(1,:,speedInd,cohInds(i)),...
            transitionInd,timeconstant,[1,1]);
        input(2,:,i) = inputStreamInterpolate(inputInitCoh(2,:,speedInd,cohInds(1)),...
            inputInitCoh(2,:,speedInd,cohInds(i)),...
            transitionInd,timeconstant,inputGain);
    end
    % input(1,:,:) = inputInitCoh(1,:,speedInd,cohInds);
    input(3,:,:) = inputInitCoh(3,:,3,1:3);
end

%% Run simulation of DynCoh
for i = 1:length(cohInds)
    [~,kappas(:,:,i)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
        't',modelFEF.t,'us',input(:,:,i),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
end


%% Select data from neuronTyping related to this subject


switch dynDataLoad
    case 'neuronTyping'
        neuralData = load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240531.mat','Rdyn','cellID','subjects');
        subjecti = find(strcmp(neuralData.subjects,subject));
    case 'dcpObjects'
        switch subject
            case 'ar'
                temp = load('~/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat');
                dcp = temp.dcp;
                sourceDirectory = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle';
                
            case 'fr'
                temp = load('~/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat');
                dcp = temp.dcp;
                sourceDirectory = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick';
                
        end
        
        [Rdyn,cellID,passCutoff,N] = calcRdynLocally(dcp,sourceDirectory,cohLookUpTable,combo);
        
        Rdyn = Rdyn(:,:,passCutoff);
        cellID = cellID(passCutoff,:,:);
        N = N(passCutoff,:);
        
%         m = squeeze(max(neuronTyping.Rinit,[],[1,2,3]))*1000;
%         m2 = squeeze(max(Rdyn,[],[1,2]))*1000;
%         Rdyn = Rdyn(:,:,m<=150 & m2<=150);
%         cellID = cellID(m<=150 & m2<=150,:,:);
%         N = N(m<=150 & m2<=150,:);
        neuralData.Rdyn = Rdyn;
        subjecti = find(strcmp(neuronTyping.subjects,subject));
        cellID(:,:,4) = subjecti;
        neuralData.cellID = cellID;
        neuralData.N = N;
end

[~,ia,ib] = intersect(squeeze(neuronTyping.cellID(:,1,:)),squeeze(neuralData.cellID(:,1,:)),'rows');
temp = false(size(neuronTyping.cellID,1),1);
temp(ia) = true;
LIA = temp & neuronTyping.cellID(:,1,4) == find(strcmp(neuronTyping.subjects,subject));
data_t = neuronTyping.initCoh.neuron_t;
Y = neuronTyping.Y(LIA & neuronTyping.cellID(:,1,4),:);
idx = neuronTyping.idx(LIA);

[~,ia] = ismember(squeeze(neuronTyping.cellID(neuronTyping.cellID(:,1,4)==find(strcmp(neuronTyping.subjects,subject)),1,:)),...
    squeeze(neuralData.cellID(neuralData.cellID(:,1,4) == subjecti,1,:)),'rows');
temp = false(size(neuronTyping.cellID(neuronTyping.cellID(:,1,4)==find(strcmp(neuronTyping.subjects,subject)),1,:),1),1);
temp(find(ia)) = true;
LIA2 = temp;

Rinit = neuronTyping.Rinit(:,:,:,neuronTyping.cellID(:,1,4) == subjecti);
Rinit = Rinit(:,:,:,LIA2);

temp = false(size(neuralData.cellID,1),1);
temp(ib) = true;
LIB = temp & neuralData.cellID(:,1,4) == subjecti;
R = neuralData.Rdyn(:,:,LIB);

if isfield(neuralData,'N') && trialCutoff > 0
    N = neuralData.N(LIB,end);
    Y = Y(N>=trialCutoff,:);
    R = R(:,:,N>=trialCutoff);
    Rinit = Rinit(:,:,:,N>=trialCutoff);
end

initNormalization = max(abs(Rinit),[],[1,2,3]);
Rinit_c = Rinit./initNormalization;
Rinit_c = Rinit_c-mean(Rinit_c(data_t<=0,:,:,:),[1,2,3]);
R_c = R./permute(initNormalization,[2,3,4,1]);
R_c = R_c-nanmean(R_c(data_t<=0,:,:),[1,2]);

kappas2 = repmat(permute(kappas,[1,4,2,3]),[1,size(Y,1),1,1]);


%% Predict the firing rates and compare to actual firing rates
X = nan(size(input,2)*size(input,3),size(input,1)+size(kappas2,1)+1);
Xout = nan(size(X,1),size(R,3));

ind = 1;
for ci = 1:size(input,3)
    X(ind:ind+size(input,2)-1,1:size(input,1)) = permute(input(:,:,ci),[2,1,3]);
    X(ind:ind+size(input,2)-1,size(input,1)+1:end-1) = permute(kappas2(:,1,:,ci),[3,1,2,4]);
    X(ind:ind+size(input,2)-1,end) = 1;
    Xout(ind:ind+size(input,2)-1,:) = squeeze(R_c(data_t<=t(end),cohLookUpTable(ci),:));
    ind = ind + size(input,2);
end

Btemp = predictFiringRates.B(:,LIA2);
if isfield(neuralData,'N') && trialCutoff > 0
    N = neuralData.N(LIB,end);
    Btemp = Btemp(:,N>=trialCutoff);
end
Xhat = X*Btemp;

dXout = Xout - repmat(Xout(1:length(t),:),[length(cohLookUpTable),1]); %
dXhat = Xhat - repmat(Xhat(1:length(t),:),[length(cohLookUpTable),1]); %repmat(control,[size(input,3),1]); %

control = squeeze(Rinit_c(data_t<=t(end),2,2,:));

Rhat = nan([length(t),size(R_c,3),size(input,3)]);
Rout = nan(size(Rhat));
dR = nan(size(Rhat));
dRhat = nan(size(Rhat));
for ci = 1:size(input,3)

    Rout(:,:,ci) = Xout((ci-1)*length(t)+1:(ci)*length(t),:);
    Rhat(:,:,ci) = Xhat((ci-1)*length(t)+1:(ci)*length(t),:);
    dR(:,:,ci) = dXout((ci-1)*length(t)+1:(ci)*length(t),:);
    dRhat(:,:,ci) = dXhat((ci-1)*length(t)+1:(ci)*length(t),:);
end
    ci = [2,3];
    CC = corrcoef([dR(t>correlationWin(1) & t<correlationWin(2),:,ci(1)),dRhat(t>correlationWin(1) & t<correlationWin(2),:,ci(1));...
        dR(t>correlationWin(1) & t<correlationWin(2),:,ci(2)),dRhat(t>correlationWin(1) & t<correlationWin(2),:,ci(2))]);
    cc = diag(CC(size(dR,2)+1:end,1:size(dR,2)));
% end

dcontrol = Xout(1:length(t),:) - control;
for ci = 1:size(input,3)
    snr(ci,:) = mean(dR(t>averagingWin(1) & t<averagingWin(2),:,ci))./sqrt(var(dcontrol(t>averagingWin(1) & t<averagingWin(2),:)) + var(dR(t>averagingWin(1) & t<averagingWin(2),:,ci)));
end
significantResponse = (abs(snr(2,:)) > 1.75 | abs(snr(3,:)) > 1.75);


%% Plot
c = lines;
hRR = figure('Name','Reduced rank model output to coherence pulse','Position', [363 445 1697 589]);
% hInputs = figure('Name','Inputs','Position', [368 450 1687 268]);
for inputi = 1:size(input,1)
    subplot(2,size(input,1),inputi)
    for ci = 1:length(cohInds)
        plot(t,squeeze(inputInitCoh(inputi,:,speedInd,cohInds(ci))),'Color',1-cohs(cohInds(ci))*ones(1,3)/100,...
            'DisplayName',['Speed = 10 deg/s, coherence = ' num2str(cohs(cohInds(ci))) '%'])
        hold on
    end
%     ind = 1:size(inputInitCoh,3);
    for ci = 2:length(cohInds)
        plot(t,squeeze(input(inputi,:,ci)),'Color',c(ci,:),...
            'DisplayName',['Speed = 10 deg/s, \Delta Coherence = ' num2str(deltaCoh_exp(ci)) '%'])
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Input dim ' num2str(inputi)])
    ylim([-1,1])
end

for ki = 1:size(kappas,1)
    subplot(2,size(input,1),size(input,1)+ki)
    for ci = 1:size(rrModel.kappas,4)
        plot(t,rrModel.kappas(ki,:,speedInd,ci),'Color',1-cohs(ci)*ones(1,3)/100,...
            'DisplayName',['Speed = 10 deg/s, Coh = ' num2str(cohs(ci))])
        hold on
    end
    for ci = 2:length(cohInds)
        plot(t,kappas(ki,:,ci),'Color',c(ci,:),...
            'DisplayName',['Speed = 10 deg/s, \Delta Coherence = ' num2str(deltaCoh_exp(ci)) '%'])
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Activity along ' rrModel.modelFEF.dimNames{ki} ' mode'])
end


%% Prediction of neural firing response to perturbations

hMeanChanges = figure('Name','Mean changes','Position',[532 603 1044 427]);
subplot(1,2,1)
histogram(mean(dcontrol(t>averagingWin(1) & t<averagingWin(2),:)))
xlabel('Control change from initCoh data')
ylabel('# of neurons')

subplot(1,2,2)
for ci = 2:size(input,3)
    x = mean(dR(t>averagingWin(1) & t<averagingWin(2),:,ci));
    y = mean(dRhat(t>averagingWin(1) & t<averagingWin(2),:,ci));
    scatter(x(:),y(:),'filled',...
        'MarkerFaceColor',c(ci,:),'MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2,...
        'DisplayName',['\Delta Coherence = ' num2str(deltaCoh_exp(ci)) '%'])
    hold on
%     axis([-1, 1, -1, 1])
%     axis([-1, 1, -1, 1]*0.5)
    plotUnity;
    axis square
    xlabel('Change from control')
    ylabel('Predicted change based on model')
end

hCC = figure('Name','Distribution of model correlation coefficients','Position',[1960 552 560 250]);
histogram(cc,linspace(-1,1,20),...
    'DisplayName','All neurons')
hold on
histogram(cc(significantResponse),linspace(-1,1,20),...
    'DisplayName','SNR > 1.75')
plotVertical([mean(cc(~isnan(cc))) ...
    mean(cc(~isnan(cc)))+1.96*std(cc(~isnan(cc)))/sqrt(sum(~isnan(cc))) ...
    mean(cc(~isnan(cc)))-1.96*std(cc(~isnan(cc)))/sqrt(sum(~isnan(cc)))]);
lineProps.Color = c(2,:);
plotVertical([mean(cc(significantResponse)) ...
    mean(cc(significantResponse))+1.96*std(cc(significantResponse))/sqrt(sum(significantResponse)) ...
    mean(cc(significantResponse))-1.96*std(cc(significantResponse))/sqrt(sum(significantResponse))],...
    'lineProperties',lineProps);
xlabel('r value')
ylabel('# of neurons')

hExNeuron = figure('Name','Example neuron');
cellID = squeeze(neuronTyping.cellID(neuronTyping.cellID(:,1,4)==find(strcmp(neuronTyping.subjects,subject)),1,:));
cellID = cellID(LIA2,:);
if isfield(neuralData,'N') && trialCutoff > 0
    N = neuralData.N(LIB,end);
    cellID = cellID(N >= trialCutoff,:);
end
switch subject
    case 'ar'
        listIndex = find(ismember(cellID(:,1:2), [77,140], 'rows'));
    case 'fr'
        listIndex = find(ismember(cellID(:,1:2), [172,11], 'rows'));
end


subplot(2,1,1)
for ci = 2:size(input,3)
    plot(t,squeeze(R(data_t<=t(end),cohLookUpTable(ci),listIndex))*1000,'Color',c(ci,:))
    hold on
end
plot(t,squeeze(R(data_t<=t(end),5,listIndex))*1000,'Color',c(5,:))
plot(t,squeeze(Rinit(data_t<=t(end),2,2,listIndex))*1000,'Color',0.4*ones(3,1))
plotVertical(motionChangeTimes(2:3));
ax = axis;
text(ax(1)+0.1*(ax(2)-ax(1)),ax(3) + 0.1*(ax(4)-ax(3)),['r value: ' num2str(cc(listIndex))])
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

subplot(2,1,2)
for ci = 2:size(input,3)
    plot(t,squeeze(dR(:,listIndex,ci)),'Color',c(ci,:))
    hold on
    plot(t,squeeze(dRhat(:,listIndex,ci)),'--','Color',c(ci,:))
end
plotVertical(correlationWin);
ax = axis;
text(ax(1)+0.1*(ax(2)-ax(1)),ax(3) + 0.1*(ax(4)-ax(3)),['r value: ' num2str(cc(listIndex))])
xlabel('Time from motion onset (ms)')
ylabel('Coherence perturbation - control (normalized)')


%% Save figures
if saveFigures
    
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' subject ...
        '/predictDynCohRRdynamics/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    
    savefig(hRR,[saveLocation '/reducedRankModelInputs&Output.fig'])
    savefig(hMeanChanges,[saveLocation '/coherencePulseResponse_predicted_v_data.fig'])
    savefig(hCC,[saveLocation '/neuronPredictivePowerDistribution.fig'])
    savefig(hExNeuron,[saveLocation '/exampleNeuronCoherencePerturbationResponse.fig'])
    
    if useMTdata
        for i = 1:length(figureHandles)
            savefig(figureHandles{i},[saveLocation '/MTfigures_' num2str(i) '.fig'])
        end
    end
end


%% Save results
if saveResults
    close all
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' subject ...
        '/ReducedRankModel'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/predictDynCohRRdynamics' datestr(now,'yyyymmdd')],'-v7.3')
    
end

%% Functions
function input = inputStreamInterpolate(stream1,stream2,transitionInd,timeconstant,inputGain)
    %%
    t = 0:100;
    f = exp(-t/timeconstant);
    f = f/sum(f);
    input = stream1;
    streamDiff = stream2(:,transitionInd:end)-stream1(:,transitionInd:end);
    for inputi = 1:size(input,1)
        streamDiff_pos = zeros(size(streamDiff(inputi,:)));
        streamDiff_pos(streamDiff(inputi,:) > 0) = streamDiff(inputi,streamDiff(inputi,:) > 0);
        streamDiff_neg = zeros(size(streamDiff(inputi,:)));
        streamDiff_neg(streamDiff(inputi,:) < 0) = streamDiff(inputi,streamDiff(inputi,:) < 0);
        input(inputi,transitionInd:end) = input(inputi,transitionInd:end) + ...
            inputGain(1)*streamDiff_pos + inputGain(2)*streamDiff_neg;
        %     input(:,transitionInd-length(t):end) = filter(f,1,input(:,transitionInd-length(t):end),[],2);
        input(inputi,:) = lowpass(input(inputi,:),1/timeconstant,1000);
    end
    
function input = inputStreamInterpolate2(stream1,stream2,transitionInd,timeconstant)
    %%
    input = stream1;
    for inputi = 1:size(input,1)
        input(inputi,transitionInd:end) = stream2(:,transitionInd:end);
        input(inputi,:) = lowpass(input(inputi,:),1/timeconstant,1000);
    end
    
function [Rdyn,cellID,passCutoff,N] = calcRdynLocally(dcp,sourceDirectory,index,combo)
    %%
    
    % Preallocate
    Rdyn = nan(1701,5,500);
    passCutoff = nan(500,1);
    cellID = nan(500,100,3);
    N = nan(500,5);
    indx = 1;
    
    % DynCoh data
    
    for filei = 1:length(dcp)
        disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
        
        dF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        dynCoh = dF.dynCoh;
        
        if ~isempty(dynCoh.r)
            
            passCutoff(indx:indx+length(dynCoh.passCutoff)-1) = dynCoh.passCutoff;
            
            
            % Get data for each neuron
                unitInd = 1:length(dynCoh.preferredDirectionRelative);
            for uniti = unitInd
                for seqi = 1:length(index)
                    [~,condLogical] = trialSort(dynCoh,dynCoh.preferredDirectionRelative(uniti),...
                        10,NaN,NaN,combo{seqi},[0,4,6,8]);
                    
                    Rdyn(:,index(seqi),indx) = nanmean(dynCoh.r(:,condLogical,uniti),2);
                    N(indx,index(seqi)) = sum(condLogical);
                end
                
                for j = 1:length(dynCoh.unitIndex)
                    cellID(indx,j,1) = filei;
                    cellID(indx,j,2) = dynCoh.unitIndex(uniti);
                    cellID(indx,j,3) = dynCoh.unitIndex(j);
                end 
                indx = indx+1;            
            end
        end
    end
    Rdyn = Rdyn(:,:,1:indx-1);
    cellID = cellID(1:indx-1,:,:);
    passCutoff = logical(passCutoff(1:indx-1));
    N = N(1:indx-1,:);