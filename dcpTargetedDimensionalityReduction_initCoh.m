function dcpTargetedDimensionalityReduction_initCoh(varargin)
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
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'dimNames',dimNames_default)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'optimize2timePoints',false)
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
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
dimNames = Parser.Results.dimNames;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
optimize2timePoints = Parser.Results.optimize2timePoints;
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
    slips = -squeeze(meanEyeSpeed(eye_t == 700,:,:)) + speeds;
    initGain(:,:,2) = initGain(:,:,2).*speeds*0.4./log2(1.4);
    initGain(:,:,3) = initGain(:,:,3).*speeds*0.4./log2(1+0.4.*speeds./slips);
end
clear init

%% Collate neural data
    
    disp('Neural analysis loop...')
    
    %% Get mean and covariance of each unit
    [Rinit, ~, cellID, passCutoff, locations] = collateFiringRates(dcp,...
        'sourceDirectory',sourceDirectory,'directions',directions,'chanMap',chanMap,...
        'rateCutoff',rateCutoff,'checkUnitType',checkUnitType,...
        'initCohCollate',true,'dynCohCollate',false);
    
    temp = load([sourceDirectory '/' dcp{1}.datapath(end-8:end-1) 'obj/initCoh' ...
             dcp{1}.datapath(end-8:end)]);
    initCoh = temp.initCoh;
    
    %% Remove data that doesn't pass cutoff
    Rinit = Rinit(:,:,:,passCutoff);
    locations = locations(passCutoff,:);
    cellID = cellID(passCutoff,:,:);
    
    %% Remove outlier rates
    m = squeeze(max(Rinit,[],[1,2,3]))*1000;
    Rinit = Rinit(:,:,:,m<=150);
    locations = locations(m<=150,:);
    cellID = cellID(m<=150,:,:);
        
    %% Mean center
    Rinit2 = Rinit;
    mRinit2 = nanmean(Rinit2,[1,2,3]);
    Rinit2 = Rinit2 - mRinit2;
        
    %% Reshape
    sz = size(Rinit2);
    Rinit3 = reshape(Rinit2,[prod(sz(1:3)),prod(sz(end))]);
            
    %% Targeted dimensionality reduction
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(Rinit3);
    zInit = Rinit2./repmat(std(Rinit2,[],[1,2,3]),[size(Rinit2,1),size(Rinit2,2),size(Rinit2,3)]);
    
    % Linear model
    Binit = nan(length(dimNames),size(zInit,4));
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
    tempGain = initGain(:,:,2);
    designMatrix2 = NaN(numel(Spds),length(dimNames));
    for di = 1:length(dimNames)
        switch dimNames{di}
            case 'Speed'
                designMatrix2(:,di) = Spds(:);
            case 'Coherence'
                designMatrix2(:,di) = Cohs(:);                
            case 'Gain'
                designMatrix2(:,di) = tempGain(:); 
            case 'Offset'
                designMatrix2(:,di) = ones(numel(Spds),1);
            otherwise
                error(['dimName ' dimNames{di} ' not recognized'])
        end
    end
    for uniti = 1:size(zInit,4)
        ztemp = reshape(zInit(initCoh.neuron_t==750,:,:,uniti),[numel(Spds),1]);
        ztemp2 = reshape(zInit(initCoh.neuron_t==150,:,:,uniti),[numel(Spds),1]);
        if optimize2timePoints
            Binit(:,uniti) = regress([ztemp; ztemp2],[designMatrix; designMatrix2]);
        else
            Binit(:,uniti) = regress(ztemp,designMatrix);
        end
    end
    
    Dtemp = nan(size(COEFF,1),size(COEFF,1));
    BinitPCA = nan(length(dimNames),size(Binit,2));
    for n = 1:24
        Dtemp(:,:,n) = COEFF(:,n)*COEFF(:,n)';
    end
    D = sum(Dtemp,3);
    BinitPCA = permute(D*Binit',[2,1,3]);
    [BinitOrth,~] = qr(BinitPCA');
    sz = size(Rinit2);
    Xinit = COEFF(:,1:10)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])';
    pInit = reshape((BinitOrth(:,1:length(dimNames))'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),length(dimNames)]);
        
        
    %% Find targeted dimensions that optimize the representation across time
    
    % Optimized decoder across time
%     taccept = initCoh.neuron_t >= 150 & initCoh.neuron_t <= 1200;
%     bOpt = nan(4,size(zInit,4));
%     bOptIsolated = nan(2,size(zInit,4));
%     for uniti = 1:size(zInit,4)
%         temp = reshape(permute(zInit(taccept,:,:,uniti),[2,3,1]),[numel(Spds)*sum(taccept),1]);
%         bOpt(:,uniti) = regress(temp,...
%             repmat([Spds(:),Cohs(:),tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
%         bOptIsolated(:,uniti) = regress(temp,...
%             repmat([tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
%     end    
%     
%     % Find unit vector 
%     normBopt = nan(size(bOpt));
%     normBoptIsolated = nan(size(bOptIsolated));
%     for dimi = 1:size(bOpt,1)
%         normBopt(dimi,:) = bOpt(dimi,:)/norm(bOpt(dimi,:));
%     end
%     normBoptIsolated(1,:) = bOptIsolated(1,:)/norm(bOptIsolated(1,:));
%     normBoptIsolated(2,:) = bOptIsolated(2,:)/norm(bOptIsolated(2,:));
%     
%     % Project into major PC dimensions
%     bOptPCA = permute(D*bOpt',[2,1]);
%     
%     % Find orthonormal basis
%     [bOptOrth,~] = qr(bOpt');
%     
%     % Project along optimal subspace dimensions
%     sz = size(zInit);
%     pBopt = reshape((bOptOrth(:,1:4)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
%     pBoptIsolated = reshape((normBoptIsolated*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),2]);
     

%% Save results file
if saveResults
    
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' dcp{1}.sname ...
        '/targetedDimensionality'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/targetedDimsInitCoh' datestr(now,'yyyymmdd')],'-v7.3')
    
end

%%
if plotOpts.On
    %% Simple avaraging
    a = Rinit(initCoh.neuron_t<=1350,:,:,:);        
    
    hGrandAverage = figure('Name','Grand average','Position',[486 733 522 420]);
    minMax = [Inf,0];
    for speedi = 1:length(speeds)
        subplot(length(speeds),1,speedi)
        for cohi = 1:length(cohs)
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speedi,cohi,:)),2)*1000,...
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
        subplot(length(speeds),1,speedi)
        ylim(minMax)
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
    end
        
    gvmrh = figure('Name',['Grand average']);
    subplot(1,2,1)
    for si = 1:length(speeds)
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
        initRatesTemp = squeeze(nanmean(a(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,:),[1,4]))*1000;
        plot(eyeSpeedTemp(si,:),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
        hold on
    end
    for si = 1:length(speeds)
        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
        initRatesTemp = squeeze(nanmean(a(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,:),[1,4]))*1000;
        plot(eyeSpeedTemp(si,:),initRatesTemp,...
            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
    end
    xlabel('Eye speed (deg/s)')
    ylabel('Spikes/s')
    axis square
    set(gca,'TickDir','out')
    
    subplot(1,2,2)
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,:),[1,4]))*1000;
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
        hold on
    end
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,:),[1,4]))*1000;
        plot(squeeze(initGain(si,:,2)),initRatesTemp,...
            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si))])
    end 
    xlabel('Behavioral gain (unitless)')
    ylabel('Spikes/s')
    axis square
    set(gca,'TickDir','out')
    
    
    %% Plot targeted dimensionality reduction
    hTargetedDimensions = figure('Name','Targeted dimension activity','Position',[63 169 522 1079]);
    
    for di = 1:length(dimNames)
        subplot(length(dimNames),1,di)
        for speedi = 1:3
            for cohi = 1:3
                plot(initCoh.neuron_t,pInit(:,speedi,cohi,di),'-','Color',[speedColors(speedi,:) cohs(cohi)/100],...
                    'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel([dimNames{di} ' related activity (a.u.)'])
    end
        
    
    %% Targeted dimension subspace
    figure('Name','Targeted subspace')
    initTaccept = initCoh.neuron_t>=-100 & initCoh.neuron_t<=1350;
    dimSelection = [2,3];
    for speedi = 1:3
        for cohi = 1:3
            plot(...
                pInit(initTaccept,speedi,cohi,dimSelection(1)),pInit(initTaccept,speedi,cohi,dimSelection(2)),...
                '-','Color',[speedColors(speedi,:) cohs(cohi)/100],...
                'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
            hold on
%             plot3(...
%                 pInit(initCoh.neuron_t==150,speedi,cohi,dimSelection(1)),pInit(initCoh.neuron_t==150,speedi,cohi,dimSelection(2)),...
%                 initGain(speedi,cohi,2),...
%                 'd','Color',[speedColors(speedi,:)],'MarkerFaceColor',[speedColors(speedi,:)],...
%                 'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
%             plot3(...
%                 pInit(initCoh.neuron_t==750,speedi,cohi,dimSelection(1)),pInit(initCoh.neuron_t==750,speedi,cohi,dimSelection(2)),...
%                 initGain(speedi,cohi,3),...
%                 'o','Color',[speedColors(speedi,:)],'MarkerFaceColor',[speedColors(speedi,:)],...
%                 'DisplayName',['Speed = ' num2str(speeds(speedi)) ', Coh = ' num2str(cohs(cohi))])
        end
    end
    grid on
    xlabel([dimNames{dimSelection(1)} ' related activity'])
    ylabel([dimNames{dimSelection(2)} ' related activity'])
    
        
    
    %% Targeted dimension activity vs gain
    gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 1000 1263]);
    for targetedDim = 1:length(dimNames)
        tempData = [];
        subplot(length(dimNames),2,2*targetedDim-1)        
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
            initRatesTemp = squeeze(nanmean(pInit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,targetedDim),1));
            plot(eyeSpeedTemp(si,:),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speeds(si))])
            hold on
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
            initRatesTemp = squeeze(nanmean(pInit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,targetedDim),1));
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
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speeds(si))])
            hold on
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
        end
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,targetedDim),1));
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
    
    
    %% Plot activity along orthonormal subspace defined by targeted optimized across time points
%     figure('Name','Activity along targeted dimensions found by optimizing over time',...
%         'Position',[73 164 522 1072])
%     for di = 1:3
%         subplot(3,1,di);
%         for si = 1:length(speeds)
%             for ci = 1:length(cohs)
%                 plot(initCoh.neuron_t,pBopt(:,si,ci,di),'-','Color',[speedColors(si,:) cohs(ci)/100])
%                 hold on
%             end
%         end
%         axis tight
%         xlabel('Time from motion onset (ms)')
%         ylabel('Activity along targeted dimension')
%         title(dimNames{di})
%         
%     end
    
    %% Plot behavioral gain vs activity along targeted dimensions
%     gvthOpt = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
%     for targetedDim = 1:3
%         tempData = [];
%         subplot(3,1,targetedDim)
%         
%         for si = 1:length(speeds)
%             initRatesTemp = squeeze(nanmean(pBopt(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,targetedDim),1));
%             plot(squeeze(initGain(si,:,3)),initRatesTemp,...
%                 'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
%                 'DisplayName',['Speed = ' num2str(speeds(si))])
%             hold on
%             tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%         end
%         for si = 1:length(speeds)
%             initRatesTemp = squeeze(nanmean(pBopt(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,targetedDim),1));
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
        
        if neuronTypingComparison.On
            savefig(hNeuronTyping,[saveLocation '/neuronTypingComparison.fig'])
        end
    end

end
