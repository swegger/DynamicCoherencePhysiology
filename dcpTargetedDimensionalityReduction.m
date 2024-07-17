function dcpTargetedDimensionalityReduction(varargin)
%%
%
%
%
%
%%

%% Defaults
plotOpts_default.On = false;
neuronTypingComparison_default.On = false;
dimNames_default = {'Speed','Coherence','Gain','Offset'};

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'pertWin',250)
addParameter(Parser,'dcpInitCohPertFile','dcpObjectsPertTemp.mat')
addParameter(Parser,'initCohPertFile',[])
addParameter(Parser,'applyLogrithmicEstimatorCorrection',false)
addParameter(Parser,'dcpDynCohFile',[])
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'sequences',[1; 2; 3; 4; 5])
addParameter(Parser,'dimNames',dimNames_default)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'neuronTypingComparison',neuronTypingComparison_default)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
directions = Parser.Results.directions;
chanMap = Parser.Results.chanMap;
pertWin = Parser.Results.pertWin;
dcpInitCohPertFile = Parser.Results.dcpInitCohPertFile;
initCohPertFile = Parser.Results.initCohPertFile;
applyLogrithmicEstimatorCorrection = Parser.Results.applyLogrithmicEstimatorCorrection;
dcpDynCohFile = Parser.Results.dcpDynCohFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
sequences = Parser.Results.sequences;
dimNames = Parser.Results.dimNames;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
neuronTypingComparison = Parser.Results.neuronTypingComparison;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;
plotOpts = Parser.Results.plotOpts;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Load dcp object
load(dcpObjectFile);

%% Find the gain from behavioral data

% DynCoh gains
disp('Determining gain from dynamic coherence trials...')
if isempty(dcpDynCohFile)
    [dyn, gain] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcp,'directions',[0,180],'pertWin',pertWin);
else
    dcpDynCoh = load(dcpDynCohFile);
    [dyn, gain] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcpDynCoh.dcp,'directions',[0,180],'pertWin',pertWin);
end
if applyLogrithmicEstimatorCorrection
    for seqi = 1:length(sequences)
        for pi = 2:length(dyn.eye.pert.t(seqi,:))
            slips(seqi,pi) = -dyn.eye.pert.m(dyn.t == dyn.eye.pert.t(seqi,pi)+100,seqi,1) + 10;
        end
    end
    gain = gain.*10*0.4./log2(1 + 0.4*10./slips);
end
eye_tDyn = dyn.t;
meanEyeSpeedDyn = dyn.eye.mean;
pertTimes = dyn.eye.pert.t;
clear dyn

% InitCoh gains
if ~isempty(initCohPertFile) && isfile(initCohPertFile)
    load(initCohPertFile,'init')
else
    disp('Determining gain from initiate coherence trials...')
    dcpInitCohPert = load(dcpInitCohPertFile);
    [init,~] = initialCohPertBehavioralAnalysis(dcp{1}.sname,'dcp',dcpInitCohPert.dcp(1:end),...
        'outlierReject',false,'win',[150 200],'keep_pert3always0deg',false,'directions',[0,180],'pertWin',pertWin);
end
eye_t = init.t;
meanEyeSpeed = init.eye.mean;
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));

if applyLogrithmicEstimatorCorrection   
    slips = -squeeze(init.eye.mean(init.t == 700,:,:)) + speeds;
    initGain(:,:,2) = initGain(:,:,2).*speeds*0.4./log2(1.4);
    initGain(:,:,3) = initGain(:,:,3).*speeds*0.4./log2(1+0.4.*speeds./slips);
end
clear init

%% Collate neural data
    
    disp('Neural analysis loop...')
    
    %% Get mean and covariance of each unit
    [Rinit, Rdyn, cellID, passCutoff, locations] = collateFiringRates(dcp,...
        'sourceDirectory',sourceDirectory,'directions',directions,'chanMap',chanMap,...
        'rateCutoff',rateCutoff,'checkUnitType',checkUnitType,...
        'initCohCollate',true,'dynCohCollate',true);
    
    temp = load([sourceDirectory '/' dcp{1}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{1}.datapath(end-8:end)]);
    initCoh = temp.initCoh;
    
    dynCohEmpty = true;
    ind = 1;
    while dynCohEmpty
        temp = load([sourceDirectory '/' dcp{ind}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{ind}.datapath(end-8:end)]);
        dynCoh = temp.dynCoh;
        if isempty(dynCoh.sequences)
            ind = ind+1;
        else
            dynCohEmpty = false;
        end
    end
    
    %% Remove data that doesn't pass cutoff
    Rinit = Rinit(:,:,:,passCutoff);
    Rdyn = Rdyn(:,:,passCutoff);
    locations = locations(passCutoff,:);
    cellID = cellID(passCutoff,:,:);
    
    %% Remove outlier rates
    m = squeeze(max(Rinit,[],[1,2,3]))*1000;
    m2 = squeeze(max(Rdyn,[],[1,2]))*1000;
    Rinit = Rinit(:,:,:,m<=150 & m2<=150);
    Rdyn = Rdyn(:,:,m<=150 & m2<=150);
    locations = locations(m<=150 & m2<=150,:);
    cellID = cellID(m<=150 & m2<150,:,:);
    
    %% Remove tail
    Rdyn = Rdyn(dynCoh.neuron_t<=1350,:,:);
    
    %% Mean center
    Rinit2 = Rinit;
    mRinit2 = nanmean(Rinit2,[1,2,3]);
    Rinit2 = Rinit2 - mRinit2;
    
    Rdyn2 = Rdyn;
    mRdyn2 = nanmean(Rdyn2,[1,2]);
    Rdyn2 = Rdyn2 - mRdyn2;
    
    %% Reshape
    sz = size(Rinit2);
    Rinit3 = reshape(Rinit2,[prod(sz(1:3)),prod(sz(end))]);
    
    sz = size(Rdyn2);
    Rdyn3 = reshape(Rdyn2,[prod(sz(1:2)),prod(sz(end))]);
        
    %% Targeted dimensionality reduction
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(Rinit3);
    zInit = Rinit2./repmat(std(Rinit2,[],[1,2,3]),[size(Rinit2,1),size(Rinit2,2),size(Rinit2,3)]);
    zDyn = Rdyn2./repmat(std(Rdyn2,[],[1,2]),[size(Rdyn2,1),size(Rdyn2,2)]);
    dynCohT = dynCoh.coh;
    dynCohT = [60*ones(sum(dynCoh.neuron_t<0),length(sequences));dynCohT];
    
    % Linear model
    Binit = nan(length(dimNames),size(zInit,4));
    Bdyn = nan(length(dimNames),size(zInit,4));
    lagT = 100;
    [Cohs, Spds] = meshgrid(cohs,speeds);
    tempGain = initGain(:,:,3);
    designMatrix = NaN(numel(Spds),length(dimNames));
    for di = 1:length(dimNames)
        switch dimNames{di}
            case 'Speed'
                designMatrix(:,di) = Spds(:);
            case 'Coherence'
                designMatrix(:,di) = Cohs(:);                
            case 'Gain'
                designMatrix(:,di) = tempGain(:); 
            case 'Offset'
                designMatrix(:,di) = ones(numel(Spds),1);
            otherwise
                error(['dimName ' dimNames{di} ' not recognized'])
        end
    end
    tiInit = find(initCoh.neuron_t==750);
    tiDyn = find(dynCoh.neuron_t == 750);
    for uniti = 1:size(zInit,4)
        ztemp = reshape(zInit(tiInit,:,:,uniti),[numel(Spds),1]);
        Binit(:,uniti) = regress(ztemp,designMatrix);
        
        ztemp = reshape(zDyn(tiDyn,:,uniti),[size(zDyn,2),1]);
        if length(unique(dynCohT(tiDyn-lagT,:))) == 1
            Bdyn(strcmp(dimNames,'Speed'),uniti) = NaN;
            Bdyn(strcmp(dimNames,'Coherence'),uniti) = NaN;
            
            if any(strcmp(dimNames,'Gain')) && any(strcmp(dimNames,'Offset'))
                temp = regress(ztemp,[gain(:,3),ones(numel(ztemp),1)]);
                Bdyn(strcmp(dimNames,'Gain'),uniti) = temp(1);
                Bdyn(strcmp(dimNames,'Offset'),uniti) = temp(2);
            elseif any(strcmp(dimNames,'Gain'))
                Bdyn(strcmp(dimNames,'Gain'),uniti) = regress(ztemp,gain(:,3));
                Bdyn(strcmp(dimNames,'Offset'),uniti) = NaN;
            elseif any(strcmp(dimNames,'Offset'))
                Bdyn(strcmp(dimNames,'Offset'),uniti) = regress(ztemp,ones(numel(ztemp),1));
                Bdyn(strcmp(dimNames,'Gain'),uniti) = NaN;
            else
                Bdyn(strcmp(dimNames,'Gain'),uniti) = NaN;
                Bdyn(strcmp(dimNames,'Offset'),uniti) = NaN;
            end
        else
            warning('off','all')
            tempDesignMatrix = nan(numel(ztemp),length(dimNames));
            designInds = 1:length(dimNames);
            for di = 1:length(dimNames)
                switch dimNames{di}
                    case 'Speed'
                        tempDesignMatrix(:,di) = NaN;
                    case 'Coherence'
                        tempDesignMatrix(:,di) = dynCohT(tiDyn-lagT,:)';
                    case 'Gain'
                        tempDesignMatrix(:,di) = gain(:,3);
                    case 'Offset'
                        tempDesignMatrix(:,di) = ones(numel(ztemp),1);
                    otherwise
                        error(['dimName ' dimNames{di} ' not recognized'])
                end
            end
            designInds = designInds(~strcmp(dimNames,'Speed'));
            temp = regress(ztemp,tempDesignMatrix(:,~strcmp(dimNames,'Speed')));
            Bdyn(designInds,uniti) = temp;
            Bdyn(strcmp(dimNames,'Speed'),uniti) = NaN;
            warning('on','all')
        end
        
    end
    
    Dtemp = nan(size(COEFF,1),size(COEFF,1));
    BinitPCA = nan(length(dimNames),size(Binit,2),size(Binit,3));
    BdynPCA = nan(length(dimNames),size(Bdyn,2),size(Bdyn,3));
    for n = 1:24
        Dtemp(:,:,n) = COEFF(:,n)*COEFF(:,n)';
    end
    D = sum(Dtemp,3);
    BinitPCA = permute(D*Binit',[2,1]);
    [BinitOrth,~] = qr(BinitPCA');
    sz = size(Rinit2);
    pInit = reshape((BinitOrth(:,1:length(dimNames))'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),length(dimNames)]);
    
    BdynPCA = permute(D*Bdyn',[2,1]);
    designInds = 1:length(dimNames);
    designInds = designInds(~strcmp(dimNames,'Speed'));
    [BdynOrth,~] = qr(BdynPCA(designInds,:)');
    sz = size(Rdyn2);
    pDynTemp = reshape((BdynOrth(:,1:sum(~strcmp(dimNames,'Speed')))'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),sum(~strcmp(dimNames,'Speed'))]);
    pDyn = nan(size(zDyn,1),size(zDyn,2),length(dimNames));
    pDyn(:,:,designInds) = pDynTemp;
    pDynCross = reshape((BinitOrth(:,1:length(dimNames))'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),length(dimNames)]);
    
        
    %% Find targeted dimensions that optimize the representation across time
%     
%     % Optimized decoder across time
%     taccept = initCoh.neuron_t >= 150 & initCoh.neuron_t <= 1200;
%     bOpt = nan(4,size(zInit,4));
%     bOptIsolated = nan(2,size(zInit,4));
%     bOptDynCoh = nan(3,size(zDyn,3));
%     for uniti = 1:size(zInit,4)
%         temp = reshape(permute(zInit(taccept,:,:,uniti),[2,3,1]),[numel(Spds)*sum(taccept),1]);
%         bOpt(:,uniti) = regress(temp,...
%             repmat([Spds(:),Cohs(:),tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
%         bOptIsolated(:,uniti) = regress(temp,...
%             repmat([tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
%     end
%     
%     taccept = taccept(dynCoh.neuron_t <=1350);
%     for uniti = 1:size(zDyn,3)
%         temp = reshape(permute(zDyn(taccept,:,uniti),[2,1]),[size(zDyn,2)*sum(taccept),1]);
%         bOptDynCoh(:,uniti) = regress(temp,...
%             [reshape(permute(dynCoh.coh(taccept,:),[2,1]),[size(zDyn,2)*sum(taccept),1]) repmat(gain(:,3),[sum(taccept),1]) ones(sum(taccept)*size(zDyn,2),1)]);
%     end
%     
%     
%     % Find unit vector 
%     normBopt = nan(size(bOpt));
%     normBoptDynCoh = nan(size(bOptDynCoh));
%     normBoptIsolated = nan(size(bOptIsolated));
%     for dimi = 1:size(bOpt,1)
%         normBopt(dimi,:) = bOpt(dimi,:)/norm(bOpt(dimi,:));
%     end
%     for dimi = 1:size(bOptDynCoh,1)
%         normBoptDynCoh(dimi,:) = bOptDynCoh(dimi,:)/norm(bOptDynCoh(dimi,:));
%     end
%     normBoptIsolated(1,:) = bOptIsolated(1,:)/norm(bOptIsolated(1,:));
%     normBoptIsolated(2,:) = bOptIsolated(2,:)/norm(bOptIsolated(2,:));
%     
%     % Project into major PC dimensions
%     bOptPCA = permute(D*bOpt',[2,1]);
%     bOptDynCohPCA = permute(D*bOptDynCoh',[2,1]);
%     
%     % Find orthonormal basis
%     [bOptOrth,~] = qr(bOpt');
%     [bOptDynCohOrth,~] = qr(bOptDynCoh');
%     
%     % Project along optimal subspace dimensions
%     sz = size(zInit);
%     pBopt = reshape((bOptOrth(:,1:4)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
%     pBoptIsolated = reshape((normBoptIsolated*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),2]);
%     
%     sz = size(zDyn);
%     pOptDynCross = reshape((bOptOrth(:,1:4)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),4]);
%     pOptDyn = reshape((bOptDynCohOrth(:,1:3)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),3]);
%     

%% Save results file
if saveResults
    
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' dcp{1}.sname ...
        '/targetedDimensionality'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/targetedDimsDynCohAndInitCoh' datestr(now,'yyyymmdd')],'-v7.3')
    
end

%%
if plotOpts.On
    %% Simple avaraging
    a = Rinit(dynCoh.neuron_t<=1350,:,:,:);
    b = nan(size(Rdyn));
    c = nan(size(Rdyn));
    for seqi = 1:size(Rdyn,2)
        t20 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 20);
        if ~isempty(t20)
            b(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:) = ...
                Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:) - ...
                Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),5,:);
        end
        
        t100 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 100);
        if ~isempty(t100)
            c(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:) = ...
                Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:) - ...
                Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),5,:);
        end
    end
        
    
    hGrandAverage = figure('Name','Grand average','Position',[486 733 1650 420]);
    minMax = [Inf,0];
    for speedi = 1:length(speeds)
        subplot(length(speeds),2,1+(speedi-1)*2)
        for cohi = 1:length(cohs)
            plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speedi,cohi,:)),2)*1000,...
                'Color',[speedColors(speedi,:), cohs(cohi)/100],...
                'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
            hold on
        end
        temp = ylim;
        if minMax(1) > temp(1)
            minMax(1) = temp(1);
        end
        if minMax(2) < temp(2)
            minMax(2) = temp(2);
        end
    end
    for speedi = 1:length(speeds)
        subplot(length(speeds),2,1+(speedi-1)*2)
        ylim(minMax)
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
    end
    
    speedi = find(speeds == 10);
    subplot(length(speeds),2,2 + (speedi-1)*2)
    for seqi = 1:length(sequences)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,seqi,:)),2)*1000,'Color',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi)])
        hold on
    end
    ylim(minMax)
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    
    gvmrh = figure('Name',['Grand average']);    
    subplot(1,2,1)
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 450 & eye_tDyn <= 450,:,:),1));
    dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,:),[1,3]))*1000;
    for seqi = 1:length(sequences)
        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi) ', Pert ' num2str(pertTimes(seqi,2)) ' ms'])
        hold on
    end
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 750 & eye_tDyn <= 750,:,:),1));
    dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,:),[1,3]))*1000;
    for seqi = 1:length(sequences)
        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi) ', Pert ' num2str(pertTimes(seqi,3)) ' ms'])
    end    
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 1050 & eye_tDyn <= 1050,:,:),1));
    dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,:),[1,3]))*1000;
    for seqi = 1:length(sequences)
        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi) ', Pert ' num2str(pertTimes(seqi,4)) ' ms'])
    end
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:),[1,4]))*1000;
        plot(eyeSpeedTemp(si,:),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
    end
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,:),[1,4]))*1000;
        plot(eyeSpeedTemp(si,:),initRatesTemp,...
            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
    end
    axis square
    set(gca,'TickDir','out') 
    xlabel('Eye speed (deg/s)')
    ylabel('Spikes/s')
    
    subplot(1,2,2)
    dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,:),[1,3]))*1000;
    for seqi = 1:length(sequences)
        plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi) ', Pert ' num2str(pertTimes(seqi,2)) ' ms'])
        hold on
    end
    dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,:),[1,3]))*1000;
    for seqi = 1:length(sequences)
        plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi) ', Pert ' num2str(pertTimes(seqi,3)) ' ms'])
    end
    dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,:),[1,3]))*1000;
    for seqi = 1:length(sequences)
        plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi) ', Pert ' num2str(pertTimes(seqi,4)) ' ms'])
    end
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:),[1,4]))*1000;
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
    end
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,:),[1,4]))*1000;
        plot(squeeze(initGain(si,:,2)),initRatesTemp,...
            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
    end
    axis square
    set(gca,'TickDir','out') 
    xlabel('Behavioral gain (unitless)')
    ylabel('Spikes/s')
    
    
    %% Plot targeted dimensionality reduction
    
    hTargetedDimensions = figure('Name','Targeted dimension activity','Position',[63 169 1606 1079]);
    
    for di = 1:length(dimNames)
        optH(1+(di-1)*3) = subplot(length(dimNames),3,1+(di-1)*3);
        for si = 1:length(speeds)
            for ci = 1:length(cohs)
                plot(initCoh.neuron_t,pInit(:,si,ci,di),'-','Color',[speedColors(si,:) cohs(ci)/100],...
                    'DisplayName',['Speed = ' num2str(speeds(si)) ', Coh = ' num2str(cohs(ci))])
                hold on
            end
        end
        axis tight
        xlabel('Time from motion onset (ms)')
        ylabel('Activity along targeted dimension')
        title(dimNames{di})
        
        optH(2+(di-1)*3) = subplot(length(dimNames),3,2+(di-1)*3);
        for seqi = 1:length(sequences)
            plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,di),'-','Color',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
            hold on
        end
        axis tight
        xlabel('Time from motion onset (ms)')
        ylabel('Activity along targeted dimension')
        title(dimNames{di})
        
        optH(3+(di-1)*3) = subplot(length(dimNames),3,3+(di-1)*3);
        for seqi = 1:length(sequences)
            plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDyn(:,seqi,di),'-','Color',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
            hold on
        end
        axis tight
        xlabel('Time from motion onset (ms)')
        ylabel('Activity along targeted dimension')
        title(dimNames{di})
    end
        
    linkaxes(optH,'xy')
    
    
    
    %% Targeted dimension subspace
    figure('Name','Targeted subspace')
    subplot(1,3,1)
    initTaccept = initCoh.neuron_t>=-100 & initCoh.neuron_t<=1350;
    dynTaccept = dynCoh.neuron_t>=-100 & dynCoh.neuron_t<=1350;
    dimSelection = [1,3];
    for speedi = 1:length(speeds)
        for cohi = 1:length(cohs)
            plot(...
                pInit(initTaccept,speedi,cohi,dimSelection(1)),pInit(initTaccept,speedi,cohi,dimSelection(2)),...
                '-','Color',[speedColors(speedi,:) cohs(cohi)/100],...
                'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
            hold on
        end
    end
    grid on
    xlabel([dimNames{dimSelection(1)} ' related activity'])
    ylabel([dimNames{dimSelection(2)} ' related activity'])
    
    subplot(1,3,2)
    for seqi = 1:5
        plot(...
            pDynCross(dynTaccept,seqi,dimSelection(1)),pDynCross(dynTaccept,seqi,dimSelection(2)),...
            '-','Color',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi)])
        hold on
    end
    grid on
    xlabel([dimNames{dimSelection(1)} ' related activity'])
    ylabel([dimNames{dimSelection(2)} ' related activity'])
    
    subplot(1,3,3)
    for seqi = 1:5
        plot(pDyn(dynTaccept,seqi,dimSelection(1)),pDyn(dynTaccept,seqi,dimSelection(2)),'-','Color',colors(seqi,:),...
            'DisplayName',['Sequence ' num2str(seqi)])
        hold on
    end
    xlabel([dimNames{dimSelection(1)} ' related activity'])
    ylabel([dimNames{dimSelection(2)} ' related activity'])
    
    %% Targeted dimension activity vs gain
    gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
    for targetedDim = 1:length(dimNames)
        tempData = [];
        subplot(length(dimNames),2,2*targetedDim-1)
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 450 & eye_tDyn <= 450,:,:),1));
        dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
            hold on
        end
        tempData = [tempData; eyeSpeedTemp(:),dynRatesTemp(:)];        
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 750 & eye_tDyn <= 750,:,:),1));
        dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
        end
        tempData = [tempData; eyeSpeedTemp(:),dynRatesTemp(:)];
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 1050 & eye_tDyn <= 1050,:,:),1));
        dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
        end
        tempData = [tempData; eyeSpeedTemp(:),dynRatesTemp(:)];
        
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,targetedDim),1));
            plot(eyeSpeedTemp(si,:),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speeds(si))])
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,targetedDim),1));
            plot(eyeSpeedTemp(si,:),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speeds(si))])
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        
        speed_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
        axis square
        ax = axis;
        text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(speed_projection_cc(1,2).^2)])
        xlabel('Eye speed (deg/s)')
        ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
        
        tempData = [];
        subplot(length(dimNames),2,2*targetedDim)
        dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
            hold on
        end
        tempData = [tempData; gain(:,2),dynRatesTemp(:)];
        dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
        end
        tempData = [tempData; gain(:,3),dynRatesTemp(:)];
        dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
        end
        tempData = [tempData; gain(:,4),dynRatesTemp(:)];
        
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speeds(si))])
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
        end
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speeds(si))])
            tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
        end
        
        gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
        axis square
        ax = axis;
        text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
        xlabel('Behavioral gain (unitless)')
        ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
    end
    
    
    gvth2 = figure('Name',['Behavioral gain vs activity on targeted dimension (dynCoh)'],'Position',[1956 59 570 1263]);
    for targetedDim = 1:length(dimNames)
        tempData = [];
        subplot(length(dimNames),1,targetedDim)
        dynRatesTemp = squeeze(nanmean(pDyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
            hold on
        end
        tempData = [tempData; gain(:,2),dynRatesTemp(:)];
        dynRatesTemp = squeeze(nanmean(pDyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
        end
        tempData = [tempData; gain(:,3),dynRatesTemp(:)];
        dynRatesTemp = squeeze(nanmean(pDyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
        for seqi = 1:length(sequences)
            plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                'DisplayName',['Sequence ' num2str(seqi)])
        end
        tempData = [tempData; gain(:,4),dynRatesTemp(:)];
        
        gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
        axis square
        ax = axis;
        text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
        xlabel('Behavioral gain (unitless)')
        ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
    end
    
    %% Plot activity along orthonormal subspace defined by targeted optimized across time points
%     figure('Name','Activity along targeted dimensions found by optimizing over time',...
%         'Position',[73 164 1616 1072])
%     dimensionLabels = {'Speed','Coherence','Gain','Offset'};
%     for di = 1:4
%         optH(1+(di-1)*3) = subplot(4,3,1+(di-1)*3);
%         for si = 1:length(speeds)
%             for ci = 1:length(cohs)
%                 plot(initCoh.neuron_t,pBopt(:,si,ci,di),'-','Color',[speedColors(si,:) cohs(ci)/100],...
%                     'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
%                 hold on
%             end
%         end
%         axis tight
%         xlabel('Time from motion onset (ms)')
%         ylabel('Activity along targeted dimension')
%         title(dimensionLabels{di})
%         
%         optH(2+(di-1)*3) = subplot(4,3,2+(di-1)*3);
%         for seqi = 1:length(sequences)
%             plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDynCross(:,seqi,di),'-','Color',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%             hold on
%         end
%         axis tight
%         xlabel('Time from motion onset (ms)')
%         ylabel('Activity along targeted dimension')
%         title(dimensionLabels{di})
%     end
%     
%     optH(6) = subplot(4,3,6);
%     for seqi = 1:length(sequences)
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDyn(:,seqi,1),'-','Color',colors(seqi,:),...
%             'DisplayName',['Sequence ' num2str(seqi)])
%         hold on
%     end
%     axis tight
%     xlabel('Time from motion onset (ms)')
%     ylabel('Activity along targeted dimension')
%     title(dimensionLabels{2})
%     
%     optH(9) = subplot(4,3,9);
%     for seqi = 1:length(sequences)
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDyn(:,seqi,2),'-','Color',colors(seqi,:),...
%             'DisplayName',['Sequence ' num2str(seqi)])
%         hold on
%     end
%     axis tight
%     xlabel('Time from motion onset (ms)')
%     ylabel('Activity along targeted dimension')
%     title(dimensionLabels{3})
%     
%     optH(12) = subplot(4,3,12);
%     for seqi = 1:length(sequences)
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDyn(:,seqi,3),'-','Color',colors(seqi,:),...
%             'DisplayName',['Sequence ' num2str(seqi)])
%         hold on
%     end
%     axis tight
%     xlabel('Time from motion onset (ms)')
%     ylabel('Activity along targeted dimension')
%     title(dimensionLabels{3})
%     
%     linkaxes(optH,'xy')
%     
    %% Plot behavioral gain vs activity along targeted dimensions
%     gvthOpt = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
%     for targetedDim = 1:4
%         tempData = [];
%         subplot(4,1,targetedDim)
%         dynRatesTemp = squeeze(nanmean(pOptDynCross(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
%         for seqi = 1:length(sequences)
%             plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%             hold on
%         end
%         tempData = [tempData; gain(:,2),dynRatesTemp(:)];
%         dynRatesTemp = squeeze(nanmean(pOptDynCross(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
%         for seqi = 1:length(sequences)
%             plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%         end
%         tempData = [tempData; gain(:,3),dynRatesTemp(:)];
%         dynRatesTemp = squeeze(nanmean(pOptDynCross(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
%         for seqi = 1:length(sequences)
%             plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%         end
%         
%         for si = 1:length(speeds)
%             initRatesTemp = squeeze(nanmean(pBopt(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,targetedDim),1));
%             plot(squeeze(initGain(si,:,3)),initRatesTemp,...
%                 'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
%                 'DisplayName',['Speed = ' num2str(speeds(si))])
%             tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%         end
%         for si = 1:length(speeds)
%             initRatesTemp = squeeze(nanmean(pBopt(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,targetedDim),1));
%             plot(squeeze(initGain(si,:,2)),initRatesTemp,...
%                 'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
%                 'DisplayName',['Speed = ' num2str(speeds(si))])
%             tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
%         end
%             
%         gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
%         axis square
%         ax = axis;
%         text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
%         xlabel('Behavioral gain (unitless)')
%         ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
%     end
%     
%     
%     gvthOpt2 = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (dynCoh)'],'Position',[1956 59 570 1263]);
%     for targetedDim = 2:4
%         tempData = [];
%         subplot(4,1,targetedDim)
%         dynRatesTemp = squeeze(nanmean(pOptDyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim-1),1));
%         for seqi = 1:length(sequences)
%             plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%             hold on
%         end
%         tempData = [tempData; gain(:,2),dynRatesTemp(:)];
%         dynRatesTemp = squeeze(nanmean(pOptDyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim-1),1));
%         for seqi = 1:length(sequences)
%             plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%         end
%         tempData = [tempData; gain(:,3),dynRatesTemp(:)];
%         dynRatesTemp = squeeze(nanmean(pOptDyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim-1),1));
%         for seqi = 1:length(sequences)
%             plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
%                 'DisplayName',['Sequence ' num2str(seqi)])
%         end
%         tempData = [tempData; gain(:,4),dynRatesTemp(:)];
%         
%         gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
%         axis square
%         ax = axis;
%         text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
%         xlabel('Behavioral gain (unitless)')
%         ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
%     end
%     
    %% Neuron typing comparison
    if neuronTypingComparison.On
        neuronTyping = load(neuronTypingComparison.file,'cellID','Y','subjects');
        
        subjectInd = find(strcmp(neuronTyping.subjects,dcp{1}.sname));
        
        Y = neuronTyping.Y;
        Y = Y(neuronTyping.cellID(:,1,4) == subjectInd,:);
        
        hNeuronTyping = figure('Name','Comparison to neuron typing','Position',[453 902 2068 420]);
        cax = [Inf, -Inf];
        for dimi = 1:length(dimNames)
            subplot(1,length(dimNames),dimi)
            scatter(Y(:,1),Y(:,2),abs(BinitOrth(:,dimi))*1000,BinitOrth(:,dimi),'filled')
            caxTemp = caxis;
            cax(1) = min([cax(1) caxTemp(1)]);
            cax(2) = max([cax(2) caxTemp(2)]);
        end
        for dimi = 1:length(dimNames)
            subplot(1,length(dimNames),dimi)
            axis square
            caxis(cax)
            xlabel('tSNE_1')
            ylabel('tSNE_2')
            title(dimNames{dimi})
            set(gca,'TickDir','out')
        end
            
    end
    
    
    %% Save figures
    if saveFigures
        
        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' dcp{1}.sname ...
            '/targetedDimensionality/' datestr(now,'yyyymmdd')];
        if ~exist(saveLocation,'dir')
            mkdir(saveLocation)
        end
        
        savefig(hGrandAverage,[saveLocation '/grandAveragePSTH.fig'])
        savefig(gvmrh ,[saveLocation '/grandAverageVsGain.fig'])
        savefig(hTargetedDimensions ,[saveLocation '/targetedDimensionActivity.fig'])
        savefig(gvth ,[saveLocation '/targetedDimensionVsGain.fig'])
        savefig(gvth2 ,[saveLocation '/targetedDimensionVsGainDynCoh.fig'])        
        
        if neuronTypingComparison.On
            savefig(hNeuronTyping,[saveLocation '/neuronTypingComparison.fig'])
        end
        
    end

end
