function modelFEF = fitTheoreticalInputModel(dcp,varargin)
%% neuronPartialCorrelation
%
%
%%

%% Defaults
plotOpts_default.On = false;
saveOpts_default.On = false;
speedPrefOpts_default.tWin = [40,120];
speedPrefOpts_default.P0 = [16,1];
speedPrefOpts_default.ub = [254,128];
speedPrefOpts_default.lb = [0, 0];
speedPrefOpts_default.c = NaN;
speedPrefOpts_default.s = NaN;
speedPrefOpts_default.d = 0;
trainCondition_default = [true, false, true;
                          true, false, true;
                          true, false, true];
theoretical_default.weightTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';
compToBehavioralGain_default.On = false;
equalizeInputsPriorToStimulusOnset_default.On = false;
f_default = @(x)(x);

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'trainCondition',trainCondition_default)
addParameter(Parser,'P0',NaN)
addParameter(Parser,'ub',NaN)
addParameter(Parser,'lb',NaN)
addParameter(Parser,'f',f_default)
addParameter(Parser,'tWin',[-100 900])
addParameter(Parser,'tau',20)
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'weightPrior',NaN)
addParameter(Parser,'tauLB',NaN)
addParameter(Parser,'tauUB',NaN)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'compToBehavioralGain',compToBehavioralGain_default)
addParameter(Parser,'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset_default)
addParameter(Parser,'plotOpts',plotOpts_default)
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
directions = Parser.Results.directions;
trainCondition = Parser.Results.trainCondition;
P0 = Parser.Results.P0;
ub = Parser.Results.ub;
lb = Parser.Results.lb;
f = Parser.Results.f;
tWin = Parser.Results.tWin;
tau = Parser.Results.tau;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
weightPrior = Parser.Results.weightPrior;
tauLB = Parser.Results.tauLB;
tauUB = Parser.Results.tauUB;
theoretical = Parser.Results.theoretical;
compToBehavioralGain = Parser.Results.compToBehavioralGain;
equalizeInputsPriorToStimulusOnset = Parser.Results.equalizeInputsPriorToStimulusOnset;
plotOpts = Parser.Results.plotOpts;
saveOpts = Parser.Results.saveOpts;

if isnan(tauLB)
    tauLB = 1;
end
if isnan(tauUB)
    tauUB = Inf;
end
%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Get firing rate data
passCutoff = nan(1000,1);
Rinit = nan(1701,3,3,1000);
cellID = nan(1000,100,3);
indx = 1;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
        
    % InitCoh data
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    
    if ~isempty(initCoh.R)
        
        passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;
        
        % Get data for each neuron
        for uniti = 1:length(initCoh.preferredDirectionRelative)
            ind = find(directions == initCoh.preferredDirectionRelative(uniti));
            Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
            
            
            for j = 1:length(initCoh.unitIndex)
                cellID(indx,j,1) = filei;
                cellID(indx,j,2) = initCoh.unitIndex(uniti);
                cellID(indx,j,3) = initCoh.unitIndex(j);
            end
                        
            indx = indx+1;
        end
    end
end


Rinit = Rinit(:,:,:,1:indx-1);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);

% taccept = initCoh.neuron_t >= tWin(1) & initCoh.neuron_t <= tWin(2);

%Rinit = Rinit(taccept,:,:,:);

fef_t = initCoh.neuron_t;%(taccept);

%% Curate data
if compToBehavioralGain.On
    temp = load(compToBehavioralGain.file,'cellID');
    cellID2 = squeeze(cellID(:,1,1:2));
    cellIDcomp = squeeze(temp.cellID(:,1,1:2));
    listIndex = find(ismember(cellID2, cellIDcomp, 'rows'));   
    Rinit = Rinit(:,:,:,listIndex);
    cellID = cellID(listIndex,:,:);
else
    % Remove data that doesn't pass cutoff
    Rinit = Rinit(:,:,:,passCutoff);
    cellID = cellID(passCutoff,:,:);
    
    % Remove outlier rates
    m = squeeze(max(Rinit,[],[1,2,3]))*1000;
    Rinit = Rinit(:,:,:,m<=150);
    cellID = cellID(m<=150,:,:);
end

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
[MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
    'speedsMT',speeds,'cohsMT',cohs,'directionsMT',directionsMT,...
    'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
    'speedPrefOpts',speedPrefOpts);

mtNeuron_t = mtResults.mt{filei}.neuron_t;


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);

inputs = permute(interpolatedR,[4,1,2,3]);
if equalizeInputsPriorToStimulusOnset.On
    switch equalizeInputsPriorToStimulusOnset.method
        case 'conditions'
            inputs(:,mtNeuron_t <= equalizeInputsPriorToStimulusOnset.latency,:,:) = ...
                repmat(mean(inputs(:,mtNeuron_t <= equalizeInputsPriorToStimulusOnset.latency,:,:),[3,4]),[1,1,size(inputs,3),size(inputs,4)]);
        case 'time&conditions'
            inputs(:,mtNeuron_t <= equalizeInputsPriorToStimulusOnset.latency,:,:) = ...
                repmat(mean(inputs(:,mtNeuron_t <= equalizeInputsPriorToStimulusOnset.latency,:,:),[2,3,4]),[1,sum(mtNeuron_t<=equalizeInputsPriorToStimulusOnset.latency),size(inputs,3),size(inputs,4)]);
        otherwise
            error('Equalization method not recognized.')
    end
end
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);
spref2 = spref(~isnan(interpolatedR(1,1,1,:)));%spref(prod(any(isnan(inputs),2),[3,4])==0);
swidth2 = swidth(~isnan(interpolatedR(1,1,1,:)));%swidth(prod(any(isnan(inputs),2),[3,4])==0);

% Find mean transient to sustained ratio
t_2_s = mean(inputs(:,mtNeuron_t==60,:,3),[3,4])./mean(inputs(:,mtNeuron_t==300,:,3),[3,4]);

%% Iterate across neurons and fit data

% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,mtNeuron_t);
R = Rinit(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);
dt = mtNeuron_t(2)-mtNeuron_t(1);
t = mtNeuron_t(iMT);
inputs_z = (inputs-mean(inputs(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(inputs(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);
R_z = (R-mean(R,[1,2,3]))./std(R,[],[1,2,3]);
taccept = t >= tWin(1) & t <= tWin(2);
RtoFit = R(taccept,:,:,:);
RtoFit_z = (RtoFit-mean(RtoFit,[1,2,3]))./std(RtoFit,[],[1,2,3]);

%% Set MT weights
switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(spref2)'-mean(log2(speedsFEF))) ones(size(spref2'))];
        
    case 'optimal'
        % More complicated: 'optimal' decoder assuming zero correlations:
        % df/ds*I/(df/ds'*I*df/ds)
        s0 = speedsFEF(speedsFEF==10);
        sprefTemp = spref2;
        %     sprefTemp(sprefTemp < 1) = 1;
        swidthTemp = swidth2;
        %     swidthTemp(swidthTemp > 10) = 10;
        df = -log2(s0./sprefTemp)./(s0.*swidthTemp*log(2)).*exp(log2(s0./sprefTemp).^2./(2*swidthTemp.^2));
        uOpt = df*inv(eye(length(sprefTemp)))/(df*inv(eye(length(sprefTemp)))*df');
        uOpt(uOpt<-0.05) = min(uOpt(uOpt>-0.05));
        uOptNorm = uOpt/norm(uOpt);
        
        Atheory = [uOpt', ones(size(uOpt'))];
end

% Normalize 
Atheory = Atheory./vecnorm(Atheory);
% Atheory = Atheory(:,1);

%% Fit dynamic model
if isa(weightPrior,'function_handle')
    W0 = weightPrior(spref2);
elseif isnan(weightPrior)
    W0 = zeros(size(spref2));
else
    W0 = weightPrior;
end

if any(isnan(P0))
%     P0 = [tauLB/dt,0,randn(1,size(inputs,1))/1000];
    P0 = [tauLB/dt,0,randn(1,3)/1000];
end
if any(isnan(lb))
    lb = [tauLB/dt, -100, -Inf, -Inf, -Inf];
end
if any(isnan(ub))
    ub = [tauUB/dt, 100, Inf, Inf, Inf];
end

mask = ones(size(trainCondition));
mask(~trainCondition) = NaN;

% OPTIONS = optimset(@fmincon);
% OPTIONS.MaxFunEvals = 3e10;
% OPTIONS.MaxIterations = 1e10;
OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10);
for neuroni = 1:size(R,4)
    tic
    disp(['Neuron ' num2str(neuroni) ' of ' num2str(size(R,4))])
%     Rtemp = R_z(taccept,:,:,neuroni);
    Rtemp = RtoFit_z(:,:,:,neuroni);
    Rtemp = Rtemp.*permute(mask,[3,1,2]);
    inputsTemp = inputs_z(:,taccept,:,:);
    inputsTemp = inputsTemp.*permute(mask,[3,4,1,2]);
    R0 = nanmean(Rtemp(1,:,:),[2,3]);
    minimizer = @(P)minimizant(P,Atheory',R0,f,inputsTemp,Rtemp);
    [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
    
    modelFEF.tau(neuroni) = P(1)*dt;
    modelFEF.R0(neuroni) = R0;
    modelFEF.baseLine(neuroni) = P(2);
    modelFEF.v(neuroni) = P(3);
    modelFEF.z(:,neuroni) = P(4:5);
    modelFEF.fval(:,neuroni) = fval;
    modelFEF.exitflag(:,neuroni) = exitflag;
    modelFEF.output(:,neuroni) = output;
    modelFEF.lambda(:,neuroni) = lambda;
    modelFEF.grad(:,neuroni) = grad;
    modelFEF.hessian(:,:,neuroni) = hessian;
    toc
end

modelFEF.dt = dt;



%% Make predictions
for neuroni = 1:size(R_z,4)
    Rhat(:,:,:,neuroni) = LNmodel(Atheory',modelFEF.v(neuroni),modelFEF.z(:,neuroni)',modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),modelFEF.tau(neuroni)/modelFEF.dt,inputs_z,f);
end

%% Rescale to z-scores across time
RhatScaled = (Rhat.*std(RtoFit,[],[1,2,3]) + mean(RtoFit,[1,2,3]) - mean(R,[1,2,3]))./std(R,[],[1,2,3]);

%% Compare fit quality

for si = 1:length(speedsFEF)
    rvalAll = corrcoef([squeeze(R(taccept,si,2,:)) squeeze(Rhat(taccept,si,2,:))]);
    rvals(:,si) = diag(rvalAll(size(R,4)+1:end,1:size(R,4)));
end
rvalsZ = 0.5*log( (1+rvals)./(1-rvals) );
[~,rvalsort] = sort(mean(rvals,2));

RhatTemp = Rhat(taccept,:,:,:);
Yhat = nan(size(RtoFit_z,1)*numel(trainCondition(:)),size(RtoFit_z,4));
Y = nan(size(RtoFit_z,1)*numel(trainCondition(:)),size(RtoFit_z,4));
conditions = nan(size(RtoFit_z,1)*numel(trainCondition(:)),3);

idx = 1;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        Y(idx:idx+size(RtoFit_z,1)-1,:) = squeeze(RtoFit_z(:,si,ci,:));
        Yhat(idx:idx+size(RtoFit_z,1)-1,:) = squeeze(Rhat(taccept,si,ci,:));
        conditions(idx:idx+size(RtoFit_z,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit_z,1),1]);
        idx = idx + size(RtoFit_z,1);
    end
end
rvalsFit = corrcoef([Y, Yhat]);
RhatCC = diag(rvalsFit(size(Y,2)+1:end,1:size(Y,2)));
RhatZ = 0.5*log( (1+RhatCC)./(1-RhatCC) );
[~,zvalSort] = sort(RhatZ);

%% Plotting
if plotOpts.On

    %% Fit quality
    figure
    for neuroni = 1:size(R,4)
        tempFEF = R_z(:,:,:,neuroni);
%         test = SimpleFEFmodel(modelFEF.W(:,neuroni)',modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),modelFEF.tau(neuroni)/modelFEF.dt,inputs_z);
        test = Rhat(:,:,:,neuroni);
        for ci = 1:length(cohsFEF)
            for si = 1:length(speedsFEF)
                subplot(length(cohsFEF),length(speedsFEF),si + (ci-1)*length(speedsFEF))
                plot(tempFEF(taccept,si,ci)*1000,test(taccept,si,ci)*1000,'o','Color',speedColors(si,:))
                hold on
            end
        end
    end
    
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(cohsFEF),length(speedsFEF),si + (ci-1)*length(speedsFEF))
        %     axis(1000*[min(RtoFit(:))-0.1*min(RtoFit(:)), max(RtoFit(:))+0.01*max(RtoFit(:)), ...
        %         min(RtoFit(:))-0.1*min(RtoFit(:)), max(RtoFit(:))+0.01*max(RtoFit(:))])
            axis square
            plotUnity;
            xlabel('FEF data (sp/s)')
            ylabel('Fit data (sp/s)')
        end
    end
    
    %% Fit quality distribution
    figure('Name','Fit quality distribution')
    plot(RhatZ(zvalSort),'bo')
    hold on
%     plotHorizontal(2.34/sqrt(size(Ytest,1)-size(Xfit,2)));       % Approximate pval of 0.01
    xlabel('FEF neuron #')
    ylabel('Fisher transformed r-values')
    
    %% Mean and mean fit values
    figure('Name','Comparison of means','Position',[391 274 1778 990])
    lims = [Inf,-Inf];
    for ci = 1:length(cohsFEF)
        subplot(2,length(cohsFEF),ci)
        for si = 1:length(speedsFEF)
            plot(t(taccept),mean(RtoFit_z(:,si,ci,:),4),'Color',speedColors(si,:))
            hold on
            plot(t(taccept),mean(Rhat(taccept,si,ci,:),4),'k--')
        end
        axis tight
        templim = ylim;
        lims(1) = min([lims(1),templim(1)]);
        lims(2) = max([lims(2),templim(2)]);
        xlabel('Time from motion onset (ms)')
        ylabel('Mean z-scored response')
    end
    
    for ci = 1:length(cohsFEF)
        subplot(2,length(cohsFEF),length(cohsFEF)+ci)
        for si = 1:length(speedsFEF)
            plot(t,mean(R_z(:,si,ci,:),4),'Color',speedColors(si,:))
            hold on
            plot(t,mean(RhatScaled(:,si,ci,:),4),'k--')
        end
        axis tight
        templim = ylim;
        lims(1) = min([lims(1),templim(1)]);
        lims(2) = max([lims(2),templim(2)]);
        xlabel('Time from motion onset (ms)')
        ylabel('Mean z-scored response')
    end
    
    for subploti = 1:2*length(cohsFEF)
        subplot(2,length(cohsFEF),subploti)
        ylim(lims)
    end
            
    
    %% Model parameters
    figure
    ah = subplot(1,3,1);
    plot(modelFEF.v(zvalSort),'o')
    ah.YDir = 'normal';
    ylabel('Neuron weight')
    xlabel('FEF neuron')
    
    subplot(1,2,2)
    plot(modelFEF.baseLine(zvalSort),'o')
    axis tight
    xlabel('FEF  neuron')
    ylabel('Baseline rate')
    
    %% Plot weigthed inputs and nonlinearity output
    figure
    lims = repmat([Inf -Inf],[size(Atheory,2)+1,1]);
    for ri = 1:size(Atheory,2)
        for ci = 1:length(cohsFEF)
            subplot(3,size(Atheory,2)+1,ci+(ri-1)*length(cohsFEF))
            for si = 1:length(speedsFEF)
                MToutputChan(:,si,ci,ri) = inputs_z(:,:,si,ci)'*Atheory(:,ri);
                plot(t,MToutputChan(:,si,ci,ri),'Color',speedColors(si,:))
                hold on
            end
            axis tight
            templims = ylim;
            lims(ri,1) = min([lims(ri,1),templims(1)]);
            lims(ri,2) = max([lims(ri,2),templims(2)]);
        end
    end
    for ci = 1:length(cohsFEF)
        subplot(3,size(Atheory,2)+1,ci+size(Atheory,2)*length(cohsFEF))
        for si = 1:length(speedsFEF)
            plot(t,f(sum(MToutputChan(:,si,ci,:),4)),'Color',speedColors(si,:))
            hold on
        end
        axis tight
        templims = ylim;
        lims(ri+1,1) = min([lims(ri+1,1),templims(1)]);
        lims(ri+1,2) = max([lims(ri+1,2),templims(2)]);
    end
    
    for ri = 1:size(Atheory,2)+1
        for ci = 1:length(cohsFEF)
            subplot(3,size(Atheory,2)+1,ci+(ri-1)*length(cohsFEF))
            ylim(lims(ri,:));
        end
    end
    
    
    
    
    %% Compare output channels to behavioral gain
    if compToBehavioralGain.On
        temp = load(compToBehavioralGain.file,'initGain');
        initGain = temp.initGain;
        
        figure('Name','MT output vs. gain')
        rankN = size(Atheory,2);
                
        for ri = 1:rankN
                subplot(1,rankN,ri)
                for si = 1:length(speedsFEF)
                    outputTemp = squeeze(nanmean(MToutputChan(t >= 700 & t <= 800,si,:,ri),1));
                    plot(squeeze(initGain(si,:,3)),outputTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                    hold on
                    outputTemp = squeeze(nanmean(MToutputChan(t >= 100 & t <= 200,si,:,ri),1));
                    plot(squeeze(initGain(si,:,2)),outputTemp,...
                        'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                end
                xlabel('Behavioral gain (unitless)')
                ylabel(['MT input along dimension ' num2str(ri)])
        end
        
    end
    
    %% Plot reconstructed firing rates projection on gain dimension
    if compToBehavioralGain.On
        temp = load(compToBehavioralGain.file,'BinitOrth','pInit');
        BinitOrth = temp.BinitOrth;
        sz = size(Rhat);
        pInit = temp.pInit(iFEF,:,:,:);
        pInitModel = reshape((BinitOrth(:,1:4)'*reshape(RhatScaled,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
        dimNames = {'Speed','Coherence','Gain','Offset'};
        
        figure('Name','Predicted response along targeted dimensions','Position',[63 169 1606 1079])
        for dimi = 1:3
            subplot(3,2,2*(dimi-1)+1)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,pInitModel(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            title('Fit predictions')
            xlabel('Time from motion onset (ms)')
            ylabel([dimNames{dimi} ' related activity (a.u.)'])
                        
            subplot(3,2,2*(dimi-1)+2)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,pInit(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            title('Data')
            xlabel('Time from motion onset (ms)')
            ylabel([dimNames{dimi} ' related activity (a.u.)'])
        end
        
        compLims = repmat([Inf,-Inf],[3,1]);
        for dimi = 1:3
            for compi = 1:2
                subplot(3,2,compi+2*(dimi-1))
                limTemp = ylim;
                compLims(dimi,1) = min([compLims(dimi,1),limTemp(1)]);
                compLims(dimi,2) = max([compLims(dimi,2),limTemp(2)]);
            end
        end
        for dimi = 1:3
            for compi = 1:2
                subplot(3,2,compi+2*(dimi-1))
                ylim(compLims(dimi,:))
            end
        end
        
        figure('Name','Systematic differencs from observed targeted dimensions','Position',[63 169 1606 1079])
        for dimi = 1:3
            subplot(3,1,dimi)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,pInit(:,speedi,cohi,dimi)-pInitModel(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            title('Difference from fit')
            xlabel('Time from motion onset (ms)')
            ylabel([dimNames{dimi} ' related activity (a.u.)'])
        end
        
        gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
        for targetedDim = 1:3
            tempData = [];
            subplot(3,1,targetedDim)
            for si = 1:length(speedsFEF)
                initRatesTemp = squeeze(nanmean(pInitModel(t >= 700 & t <= 800,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                    'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                hold on
                tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
                initRatesTemp = squeeze(nanmean(pInitModel(t >= 100 & t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
            gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
            axis square
            ax = axis;
            text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
            xlabel('Behavioral gain (unitless)')
            ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
        end
        
    end
    
    %% Relationship to Functional topography
    if compToBehavioralGain.On
        temp = load(compToBehavioralGain.file,'NumClusters','Y','idx','centInd');
        NumClusters = temp.NumClusters;
        for typei = 1:NumClusters
            typeAngle = wrapTo2Pi(atan2(temp.Y(temp.centInd == typei,2)-mean(temp.Y(:,2)),temp.Y(temp.centInd == typei,1)-mean(temp.Y(:,1))))/2+pi;
            hue = typeAngle / (2 * pi);
            saturation = ones(size(hue));
            value = ones(size(hue));
            hsv = cat(3, hue, saturation, value);
            colorWheel(typei,:) = hsv2rgb(hsv);
        end
        colorWheel(colorWheel > 1) = 1;
        colorWheel(colorWheel < 0) = 0;
        
        figure
        scatter(temp.Y(:,1),temp.Y(:,2),50,RhatZ,'filled')
%         scatter(temp.Y(:,1),temp.Y(:,2),50,q,'filled')
        colorbar
        xlabel('tSNE_1')
        ylabel('tSNE_2')
        
%         figure
%         Wtemp = W(sortInd,qsort);
%         idx = temp.idx(qsort);
%         for typei = 1:NumClusters
%             subplot(3,NumClusters,typei)
%             for si = 1:length(speedsFEF)
%                 for ci = 1:length(cohsFEF)
%                     plot(t(taccept),nanmean(RtoFit_z(:,si,ci,temp.idx==typei),4),...
%                         'Color',[speedColors(si,:) cohsFEF(ci)/100])
%                     hold on
%                 end
%             end
%             xlabel('Time from motion onset (ms)')
%             ylabel('z-score')
%             
%             subplot(3,NumClusters,NumClusters+typei)
%             imagesc(Wtemp(:,idx==typei))
%             ylabel('MT unit (sorted by speed pref)')
%             xlabel('FEF unit in cluster')
%             
%             subplot(3,NumClusters,2*NumClusters+typei)
%             Atemp = W(:,temp.idx==typei)*Vhat(temp.idx==typei,:);
%             semilogx(spref2,Atemp(:,1),'o')
%             ylabel('MT unit weight')
%             xlabel('speed preference (deg/s)')
%         end
        
    end
end


%% Saving
if saveOpts.On
    save(saveOpts.location);
end

%% Functions

%% Function to fit
function R =LNmodel(W,v,z,baseLine,R0,tau,inputs,f)
    sz = size(inputs);
    R = nan(sz(2:4));
    R(1,:,:) = repmat(R0,[1,sz(3),sz(4)]);
    for ti = 2:sz(2)
        for si = 1:sz(3)
            for ci = 1:sz(4)
                dR(ti,si,ci) = -(R(ti-1,si,ci)-baseLine) + v*f(z*W*inputs(:,ti,si,ci));
                R(ti,si,ci) = R(ti-1,si,ci) + dR(ti,si,ci)/tau;
            end
        end
    end
    
 %% Objective
function out = minimizant(P,W,R0,f,inputs,R)
    tau = P(1);
    baseLine = P(2);
    v = P(3);
    z = P(4:5);
    Rest = LNmodel(W,v,z,baseLine,R0,tau,inputs,f);
    out = nansum((R(:) - Rest(:)).^2);