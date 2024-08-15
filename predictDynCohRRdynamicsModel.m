function predictDynCohRRdynamicsModel(subject,varargin)
%%
%
%
%
%%

%% Defaults

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


%% Check for specific files as input, otherwise set default files
if any(isnan(neuronTypingFile))
    switch subject
        case 'ar'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
        case 'fr'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';            
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
neuronTyping = load(neuronTypingFile,'subjects','Y','cellID','gainRegression','idx','Rinit','NumClusters','initCoh');
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

%% Save results
if saveResults
        
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' subject ...
        '/ReducedRankModel'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/predictDynCohRRdynamics' datestr(now,'yyyymmdd')],'-v7.3')
    
end

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




%% Save figures
if saveFigures
    
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' subject ...
        '/predictDynCohRRdynamics/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    
    savefig(hRR,[saveLocation '/reducedRankModelInputs&Output.fig'])
    
    if useMTdata
        for i = 1:length(figureHandles)
            savefig(figureHandles{i},[saveLocation '/MTfigures_' num2str(i) '.fig'])
        end
    end
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