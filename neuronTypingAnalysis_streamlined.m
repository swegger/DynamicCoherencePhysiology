function neuronTypingAnalysis_streamlined(subjects,varargin)
%%
%
%
%
%
%%

%% Defaults
plotOpts_default.On = false;

embedding_default.NumDimensions = 2;
embedding_default.Perplexity = 10;

loadFunctionalTopology_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subjects')
addParameter(Parser,'dcpObjectFile',[])
addParameter(Parser,'sourceDirectory',[])
addParameter(Parser,'initCohCollate',true)
addParameter(Parser,'dynCohCollate',true)
addParameter(Parser,'embedding',embedding_default)
addParameter(Parser,'loadFunctionalTopology',loadFunctionalTopology_default)
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'chanMap',{[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24]})
addParameter(Parser,'ClusterMethod','densityClust')
addParameter(Parser,'pertWin',250)
addParameter(Parser,'dcpInitCohPertFile','dcpObjectsPertTemp.mat')
addParameter(Parser,'initSpeed',10)
addParameter(Parser,'includeEarlyInitCohPertTime',true)
addParameter(Parser,'dcpDynCohFile',{[]})
addParameter(Parser,'applyLogrithmicEstimatorCorrection',false)
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'sequences',[1; 2; 3; 4; 5])
addParameter(Parser,'NumClusters',8)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'plotOpts',plotOpts_default)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,subjects,varargin{:})

subjects = Parser.Results.subjects;
dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
initCohCollate = Parser.Results.initCohCollate;
dynCohCollate = Parser.Results.dynCohCollate;
embedding = Parser.Results.embedding;
loadFunctionalTopology = Parser.Results.loadFunctionalTopology;
directions = Parser.Results.directions;
chanMap = Parser.Results.chanMap;
ClusterMethod = Parser.Results.ClusterMethod;
pertWin = Parser.Results.pertWin;
initSpeed = Parser.Results.initSpeed;
includeEarlyInitCohPertTime = Parser.Results.includeEarlyInitCohPertTime;
dcpDynCohFile = Parser.Results.dcpDynCohFile;
applyLogrithmicEstimatorCorrection = Parser.Results.applyLogrithmicEstimatorCorrection;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
sequences = Parser.Results.sequences;
NumClusters = Parser.Results.NumClusters;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
plotOpts = Parser.Results.plotOpts;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Preallocate
passCutoff = nan(1000,1);
Rinit = nan(1701,3,3,1000);
Rdyn = nan(1701,5,1000);
locations = nan(1000,3);
cellID = nan(1000,100,4);

indx = 1;

%% Cycle through subjects
for subjecti = 1:length(subjects)
    
    %% Set up subject data
    if ~isempty(dcpObjectFile)
        temp = load(dcpObjectFile{subjecti});
    else
        switch subjects{subjecti}
            case 'ar'
                temp = load('~/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat');
            case 'fr'
                temp = load('~/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat');
            otherwise
                error('Subject not recognized!')
        end
    end
    dcp = temp.dcp;
    if ~strcmp(dcp{1}.sname,subjects{subjecti})
        error(['dcpObjectFile dcp objects are for ' dcp{1}.sname ', not ' subjects{subjecti}])
    end
    
    if ~isempty(sourceDirectory)
        sourceDirectory = sourceDirectory{subjecti};
    else
        switch subjects{subjecti}
            case 'ar'
                sourceDirectoryTemp = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle';
            case 'fr'
                sourceDirectoryTemp = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick';
            otherwise
                error('Subject not recognized!')
        end
    end
    
    %% Get mean and covariance of each unit from each object file
    [RinitTemp, RdynTemp, cellIDTemp, passCutoffTemp, locationsTemp] = collateFiringRates(dcp,...
        'sourceDirectory',sourceDirectoryTemp,'directions',directions,'chanMap',chanMap{subjecti},...
        'rateCutoff',rateCutoff,'checkUnitType',checkUnitType,...
        'initCohCollate',initCohCollate,'dynCohCollate',dynCohCollate);
    
    % Add to data matrices
    Rinit(:,:,:,indx:indx+size(RinitTemp,4)-1) = RinitTemp;
    Rdyn(:,:,indx:indx+size(RdynTemp,3)-1) = RdynTemp;
    cellID(indx:indx+size(cellIDTemp,1)-1,:,1:3) = cellIDTemp;
    cellID(indx:indx+size(cellIDTemp,1)-1,:,4) = subjecti;
    passCutoff(indx:indx+size(passCutoffTemp,1)-1,:) = passCutoffTemp;
    locations(indx:indx+size(locationsTemp,1)-1,:) = locationsTemp;
    
    indx = indx+size(locationsTemp,1);
    
    
    %% Find the gain from behavioral data    
    if dynCohCollate
        disp('Determining gain from dynamic coherence trials...')
        if isempty(dcpDynCohFile{subjecti})
            [dyn,gain{subjecti}] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcp,'directions',[0,180],'pertWin',pertWin);
        else
            dcpDynCoh = load(dcpDynCohFile{subjecti});
            [dyn,gain{subjecti}] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcpDynCoh.dcp,'directions',[0,180],'pertWin',pertWin);
        end
        eye_tDyn{subjecti} = dyn.t;
        meanEyeSpeedDyn{subjecti} = dyn.eye.mean;        
        clear dyn
    end
    
    if initCohCollate
        disp('Determining gain from initiate coherence trials...')
        switch subjects{subjecti}
            case 'ar'
                load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat', 'init')
            case 'fr'
                load('/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212.mat', 'init')
        end
        
        slips = -squeeze(init.eye.mean(init.t == 750,:,:)) + speeds;
        
        for pi = 1:3
            for ci = 1:length(cohs)
                initGain{subjecti}(:,ci,pi) = (init.eye.pert.res(:,ci,pi)-init.eye.pert.resControl(:,ci,pi))./(0.4*speeds);
                gain95CI{subjecti}(:,ci,pi) = sqrt(init.eye.pert.resSTE(:,ci,pi).^2 + init.eye.pert.resControlSTE(:,ci,pi).^2)./(0.4*speeds)*1.64;
            end
            
            if pi == 2 && applyLogrithmicEstimatorCorrection
                % Initiation
                initGain{subjecti}(:,:,pi) = initGain{subjecti}(:,:,pi).*speeds*0.4./log2(1.4);
                gain95CI{subjecti}(:,:,pi) = gain95CI{subjecti}(:,:,pi).*speeds*0.4./log2(1.4);
            elseif pi == 3 && applyLogrithmicEstimatorCorrection
                % Closed loop
                initGain{subjecti}(:,:,pi) = initGain{subjecti}(:,:,pi).*speeds*0.4./log2(1+0.4*speeds./slips);
                gain95CI(:,:,pi) = gain95CI(:,:,pi).*speeds*0.4./log2(1+0.4*speeds./slips);
            end
        end
        eye_t{subjecti} = init.t;
        meanEyeSpeed{subjecti} = init.eye.mean;
        clear init
    end
end

%%
if initCohCollate
    temp = load([sourceDirectoryTemp '/' dcp{1}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{1}.datapath(end-8:end)]);
    initCoh = temp.initCoh;
end

if dynCohCollate
    filei = 1;
    successfulDynCohLoad = false;
    while ~successfulDynCohLoad
        tempFile = [sourceDirectoryTemp '/' dcpDynCoh.dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcpDynCoh.dcp{filei}.datapath(end-8:end) '.mat'];
        if exist(tempFile,'file')
            temp = load(tempFile);
            successfulDynCohLoad = true;
        else
            filei = filei+1;
        end
    end
    dynCoh = temp.dynCoh;
end

%% 
Rinit = Rinit(:,:,:,1:indx-1);
Rdyn = Rdyn(:,:,1:indx-1);
locations = locations(1:indx-1,:);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff);
Rdyn = Rdyn(:,:,passCutoff);
locations = locations(passCutoff,:);
cellID = cellID(passCutoff,:,:);

%% Remove outlier rates
if initCohCollate
    m = squeeze(max(Rinit,[],[1,2,3]))*1000;    
else
    m = zeros(size(squeeze(max(Rinit,[],[1,2,3]))*1000));
end
if dynCohCollate
    m2 = squeeze(max(Rdyn,[],[1,2]))*1000;
else
    m2 = zeros(size(squeeze(max(Rdyn,[],[1,2]))*1000));
end
Rinit = Rinit(:,:,:,m<=150 & m2<=150);
Rdyn = Rdyn(:,:,m<=150 & m2<=150);
locations = locations(m<=150 & m2<=150,:);
cellID = cellID(m<=150 & m2<150,:,:);

%% Remove tail
if dynCohCollate
    Rdyn = Rdyn(dynCoh.neuron_t<=1350,:,:);
end

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
    
    
%% Use nonlinear dimensionality reduction to attempt to identify neuron 'types'
if loadFunctionalTopology.On
    mapTemp = load(loadFunctionalTopology.file,'Y','idx','NumClusters');
    Y = mapTemp.Y;
    idx = mapTemp.idx;
    NumClusters = mapTemp.NumClusters;
else
    if initCohCollate && dynCohCollate
        Y = tsne([Rinit3; Rdyn3]','Distance','cosine','NumDimensions',embedding.NumDimensions,'Perplexity',embedding.Perplexity);
    elseif initCohCollate
        Y = tsne(Rinit3','Distance','cosine','NumDimensions',embedding.NumDimensions,'Perplexity',embedding.Perplexity);
    elseif dynCohCollate
        Y = tsne(Rdyn3','Distance','cosine','NumDimensions',embedding.NumDimensions,'Perplexity',embedding.Perplexity);
    end
    
    switch ClusterMethod
        case 'K-means'
            NumClusters = 3;
            idx = kmeans(Y,NumClusters);
            
        case 'densityClust'
            dist = pdist2(Y,Y);
            percNeigh = 0.02;
            % 'Gauss' denotes the use of Gauss Kernel to compute rho, and
            % 'Cut-off' denotes the use of Cut-off Kernel.
            % For large-scale data sets, 'Cut-off' is preferable owing to computational efficiency,
            % otherwise, 'Gauss' is preferable in the case of small samples (especially with noises).
            kernel = 'Gauss';
            % set critical system parameters for DensityClust
            [dc, rho] = paraSet(dist, percNeigh, kernel);
            [NumClusters, idx, centInd, haloInd] = densityClust(dist, dc, rho, true);
    end
end
    
%% Distances in functional space vs physical space
for subjecti = 1:length(subjects)
    locations2 = locations(cellID(:,1,4) == subjecti,:);
    locations2(locations2(:,1)>1,:) = [locations2(locations2(:,1)>1,2), locations2(locations2(:,1)>1,1), locations2(locations2(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
    locations2(:,end) = locations2(:,end)/1000;
    Y2 = Y(cellID(:,1,4) == subjecti,:);
    
    euclidLoc{subjecti} = pdist(locations2);
    euclidY{subjecti} = pdist(Y2);
    
    % Population level distance correlation
    [locCor,locCorP] = corrcoef(euclidLoc{subjecti}(~isnan(euclidLoc{subjecti})&~isnan(euclidY{subjecti})),...
        euclidY{subjecti}(~isnan(euclidLoc{subjecti})&~isnan(euclidY{subjecti})));
    
    % K nearest neighbor distances
    Kneighbors = 10;
    nnIdx{subjecti} = knnsearch(locations2,locations2,'K',Kneighbors+1);
    for ni = 1:size(nnIdx{subjecti},1)
        mdist(ni) = mean(pdist(Y2(nnIdx{subjecti}(ni,2:end),:)));
        randIdx = [randsample(size(Y2,1),Kneighbors+1)];
        mdistRand(ni) = mean(pdist(Y2(randIdx,:)));
    end
end

thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
[~,thetaSort] = sort(thetas);

%% Perform regression
for subjecti = 1:length(subjects)
    if initCohCollate && dynCohCollate
        gainRegression(subjecti).x = nan(2*length(speeds)*length(cohs)+3*length(sequences),NumClusters);
    elseif dynCohCollate
        gainRegression(subjecti).x = nan(3*length(sequences),NumClusters);
    elseif initCohCollate
        gainRegression(subjecti).x = nan(2*length(speeds)*length(cohs),NumClusters);
    end
    for i = 1:NumClusters
        if initCohCollate
            a = Rinit(initCoh.neuron_t<=1350,:,:,idx == i);
        else
            a = nan(size(Rinit(:,:,:,idx == i)));
        end
        b = nan(size(Rdyn(:,:,idx == i)));
        c = nan(size(Rdyn(:,:,idx == i)));
        if dynCohCollate
            for seqi = 1:size(Rdyn,2)
                t20 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 20);
                if ~isempty(t20)
                    b(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:) = ...
                        Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,idx == i) - ...
                        Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),5,idx == i);
                end
                
                t100 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 100);
                if ~isempty(t100)
                    c(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:) = ...
                        Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,idx == i) - ...
                        Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),5,idx == i);
                end
            end
        end
        
        indx = 1;
        if dynCohCollate
            dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx==i),1)),2)*1000;
            gainRegression(subjecti).x(indx:indx+length(sequences)-1,i) = dynRatesTemp;
            indx = indx+length(sequences);
            dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx==i),1)),2)*1000;
            gainRegression(subjecti).x(indx:indx+length(sequences)-1,i) = dynRatesTemp;
            indx = indx+length(sequences);
            dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx==i),1)),2)*1000;
            gainRegression(subjecti).x(indx:indx+length(sequences)-1,i) = dynRatesTemp;
            indx = indx+length(sequences);
        end
        
        if initCohCollate
            for si = 1:length(speeds)
                initRatesTemp = nanmean(squeeze(nanmean(a(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,:),1)),2)*1000;
                gainRegression(subjecti).x(indx:indx+length(cohs)-1,i) = initRatesTemp;
                indx = indx+length(cohs);
            end
            if includeEarlyInitCohPertTime
                for si = 1:length(speeds)
                    initRatesTemp = nanmean(squeeze(nanmean(a(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,:),1)),2)*1000;
                    gainRegression(subjecti).x(indx:indx+length(cohs)-1,i) = initRatesTemp;
                    indx = indx+length(cohs);
                end
            end
        end
    end
    
    regressorClusters = 1:NumClusters;
    if initCohCollate && dynCohCollate
        gainRegression(subjecti).x_mean = mean(gainRegression(subjecti).x,1);
        gainRegression(subjecti).x_std = std(gainRegression(subjecti).x,[],1);
        gainRegression(subjecti).z = gainRegression(subjecti).x;
        gainRegression(subjecti).y(1:length(sequences)) = gain{subjecti}(:,2);
        gainRegression(subjecti).y(length(sequences)+1:2*length(sequences)) = gain{subjecti}(:,3);
        gainRegression(subjecti).y(2*length(sequences)+1:3*length(sequences)) = gain{subjecti}(:,4);
        for si = 1:length(speeds)
            gainRegression(subjecti).y(3*length(sequences)+(si-1)*length(cohs)+1:3*length(sequences)+si*length(cohs)) = initGain{subjecti}(si,:,3);
        end
        if includeEarlyInitCohPertTime
            for si = 1:length(speeds)
                gainRegression(subjecti).y(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+1:3*length(sequences)+3*length(cohs)+si*length(cohs)) = initGain{subjecti}(si,:,2);
            end
        end
    elseif dynCohCollate
        gainRegression(subjecti).y(1:length(sequences)) = gain{subjecti}(:,2);
        gainRegression(subjecti).y(length(sequences)+1:2*length(sequences)) = gain{subjecti}(:,3);
        gainRegression(subjecti).y(2*length(sequences)+1:3*length(sequences)) = gain{subjecti}(:,4);
        gainRegression(subjecti).z = gainRegression(subjecti).x;
    elseif initCohCollate
        for si = 1:length(speeds)
            gainRegression(subjecti).y((si-1)*length(cohs)+1:si*length(cohs)) = initGain{subjecti}(si,:,3);
        end
        if includeEarlyInitCohPertTime
            for si = 1:length(speeds)
                gainRegression(subjecti).y(3*length(cohs)+(si-1)*length(cohs)+1:3*length(cohs)+si*length(cohs)) = initGain{subjecti}(si,:,2);
            end
        end
        gainRegression(subjecti).z = gainRegression(subjecti).x;
    end
    gainRegression(subjecti).B = regress(gainRegression(subjecti).y(:),[gainRegression(subjecti).z(:,regressorClusters) ones(size(gainRegression(subjecti).x,1),1)]);
    gainRegression(subjecti).yhat = [gainRegression(subjecti).z(:,regressorClusters) ones(size(gainRegression(subjecti).x,1),1)]*gainRegression(subjecti).B;
end

%% Assign colors to clusters for plotting
for typei = 1:NumClusters
    typeAngle = wrapTo2Pi(atan2(Y(centInd == typei,2)-mean(Y(:,2)),Y(centInd == typei,1)-mean(Y(:,1))))/2+pi;
    hue = typeAngle / (2 * pi);
    saturation = ones(size(hue));
    value = ones(size(hue));
    hsv = cat(3, hue, saturation, value);
    colorWheel(typei,:) = hsv2rgb(hsv);
end
colorWheel(colorWheel > 1) = 1;
colorWheel(colorWheel < 0) = 0;

%% Save results file
if saveResults    
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' [subjects{:}] ...
        '/neuronTypingAnalysis'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/neuronTyping' datestr(now,'yyyymmdd')],'-v7.3')
    
end

%% Plotting
if plotOpts.On
    
    %% Example neurons
    exUnits = [77,140;...
        67,186;...
        69, 81;...
        70,156;...
        75, 80;...
        57,350];
    for uniti = 1:size(exUnits,1)
        
        exIndex = exUnits(uniti,:);
        cellID2 = squeeze(cellID(:,1,1:2));
        listIndex = find(ismember(cellID2, exIndex, 'rows'));
        hExUnits(uniti) = figure('Name',['Unit ' num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)}],...
            'Position',[1553 976 973 341]);
        ylims = [Inf,-Inf];
        if initCohCollate
            for speedi = 1:length(speeds)
                subplot(2,length(speeds),speedi)
                for ci = 1:length(cohs)
                    plot(initCoh.neuron_t,Rinit(:,speedi,ci,listIndex)*1000,'Color',initColors(ci,:),...
                        'DisplayName',['Neuron ' num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)} ', Speed = ' num2str(speeds(speedi)) ' deg/s, Coh = ' num2str(cohs(ci)) '%']);
                    hold on
                end
                xlabel('Time from motion onset (ms)')
                ylabel('Spikes/s')
                title('initCoh')
                axis tight
                tempLims = ylim;
                ylims(1) = min([ylims(1),tempLims(1)]);
                ylims(2) = max([ylims(2),tempLims(2)]);
            end
        end
        
        subplot(2,length(speeds),2+length(speeds))
        if dynCohCollate
            for seqi = 1:5
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,seqi,listIndex)*1000,'Color',colors(seqi,:),...
                    'DisplayName',['Neuron ' num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)} ', Sequence = ' num2str(seqi)]);
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
        title('dynCoh')
        axis tight
        tempLims = ylim;
        ylims(1) = min([ylims(1),tempLims(1)]);
        ylims(2) = max([ylims(2),tempLims(2)]);
        
        newax = [-100 1350 ylims];
        for speedi = 1:length(speeds)
            subplot(2,length(speeds),speedi)
            axis(newax)
        end
        subplot(2,length(speeds),2+length(speeds))
        axis(newax)
        text((newax(2)-newax(1))*.5+newax(1),(newax(4)-newax(3))*.1+newax(3),...
            [num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)} ])
    end
    
    %% PSTH averaged accross neurons in a cluster
    for i = 1:NumClusters
        if initCohCollate
            a = Rinit(initCoh.neuron_t<=1350,:,:,idx == i);
        else
            a = nan(size(Rinit(:,:,:,idx == i)));
        end
        b = nan(size(Rdyn(:,:,idx == i)));
        c = nan(size(Rdyn(:,:,idx == i)));
        if dynCohCollate
            for seqi = 1:size(Rdyn,2)
                t20 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 20);
                if ~isempty(t20)
                    b(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:) = ...
                        Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,idx == i) - ...
                        Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),5,idx == i);
                end
                
                t100 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 100);
                if ~isempty(t100)
                    c(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:) = ...
                        Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,idx == i) - ...
                        Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),5,idx == i);
                end
            end
        end
%         
%         figure('Name',['Cluster ' num2str(i)])
%         subplot(2,1,1)
%         for seqi = 1:5
%             plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(c(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
%             hold on
%             plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(b(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
%         end
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)-a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(1,:))
%         hold on
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)-a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
%         plotHorizontal(0);
%         plotVertical([450 750 1050]);
%         xlabel('Time from motion onset (ms)')
%         ylabel('Excess spike/s')
%         
%         subplot(2,1,2)
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
%         hold on
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,5,idx == i)),2)*1000,'Color',colors(5,:))
%         plotVertical([450 750 1050]);
%         xlabel('Time from motion onset (ms)')
%         ylabel('Spikes/s')
        
        hClusterPSTH(i) = figure('Name',['Cluster ' num2str(i)]);
        subplot(3,1,1)
        if dynCohCollate
            for seqi = 1:5
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(dynCoh.neuron_t<=1350,seqi,idx == i)),2)*1000,...
                    'Color',colors(seqi,:))
                hold on
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(dynCoh.neuron_t<=1350,seqi,idx == i)),2)*1000,...
                    'Color',colors(seqi,:))
            end
        end
        plotVertical([450 750 1050]);
        xlabel('Time from motion onset (ms)')
        ylabel('Spike/s')
        
        subplot(3,1,2)
        if initCohCollate
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)),2)*1000,'Color',initColors(1,:))
            hold on
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)),2)*1000,'Color',initColors(3,:))
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Spike/s')
        
        subplot(3,1,3)
        if initCohCollate
            for si = 1:length(speeds)
                plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,si,2,:)),2)*1000,'Color',speedColors(si,:))
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
    end
       
    %% Gain vs average response in cluster
    for i = 1:NumClusters
        gvrh(i) = figure('Name',['Cluster ' num2str(i)]);
        for subjecti = 1:length(subjects)
            subplot(2,length(subjects),subjecti)
            if dynCohCollate
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn{subjecti}(eye_tDyn{subjecti} >= 450 & eye_tDyn{subjecti} <= 450,:,:),1));
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                        hold on
                    end
                end
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn{subjecti}(eye_tDyn{subjecti} >= 750 & eye_tDyn{subjecti} <= 750,:,:),1));
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                    end
                end
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn{subjecti}(eye_tDyn{subjecti} >= 1050 & eye_tDyn{subjecti} <= 1050,:,:),1));
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                    end
                end
            end
            if initCohCollate
                for si = 1:length(speeds)
                    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed{subjecti}(eye_t{subjecti} >= 750 & eye_t{subjecti} <= 750,:,:),1));
                    initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                    plot(eyeSpeedTemp(si,:),initRatesTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 750 ms'])
                    hold on
                end
                if includeEarlyInitCohPertTime
                    for si = 1:length(speeds)
                        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed{subjecti}(eye_t{subjecti} >= 150 & eye_t{subjecti} <= 150,:,:),1));
                        initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                        plot(eyeSpeedTemp(si,:),initRatesTemp,...
                            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 50 ms'])
                    end
                end
                
            end
            xlabel('Eye speed (deg/s)')
            ylabel('Spikes/s')
            title(subjects(subjecti))
            set(gca,'TickDir','out')
            
            subplot(2,length(subjects),subjecti + length(subjects))
            if dynCohCollate
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(gain{subjecti}(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                        hold on
                    end
                end
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(gain{subjecti}(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                    end
                end
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(gain{subjecti}(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                    end
                end
            end
            if initCohCollate
                for si = 1:length(speeds)
                    initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                    plot(squeeze(initGain{subjecti}(si,:,3)),initRatesTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 750 ms'])
                    hold on
                end
                if includeEarlyInitCohPertTime
                    for si = 1:length(speeds)
                        initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                        plot(squeeze(initGain{subjecti}(si,:,2)),initRatesTemp,...
                            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 50 ms'])
                    end
                end
                
            end
            xlabel('Behavioral gain (unitless)')
            ylabel('Spikes/s')
            title(subjects(subjecti))
            set(gca,'TickDir','out')
        end
    end
    
    %% Gain decoding
    hGainDecodingFromClusters = figure('Name','Gain decoding','Position',[1633 927 888 395]);
    for subjecti = 1:length(subjects)
        subplot(length(subjects),2,1+2*(subjecti-1))
        if initCohCollate && dynCohCollate
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi),gainRegression(subjecti).yhat(seqi),...
                    'd','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+length(sequences)),gainRegression(subjecti).yhat(seqi+length(sequences)),...
                    'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+2*length(sequences)),gainRegression(subjecti).yhat(seqi+2*length(sequences)),...
                    's','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                hold on
            end
            for si = 1:length(speeds)
                for ci = 1:length(cohs)
                    plot(gainRegression(subjecti).y(3*length(sequences)+(si-1)*length(cohs)+ci),...
                        gainRegression(subjecti).yhat(3*length(sequences)+(si-1)*length(cohs)+ci),...
                        'o','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 750 ms'])
                    hold on
                end
            end
            if includeEarlyInitCohPertTime
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        plot(gainRegression(subjecti).y(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+ci),...
                            gainRegression(subjecti).yhat(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+ci),...
                            'd','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 50 ms'])
                        hold on
                    end
                end
            end
        elseif initCohCollate
            for si = 1:length(speeds)
                for ci = 1:length(cohs)
                    plot(gainRegression(subjecti).y((si-1)*length(cohs)+ci),...
                        gainRegression(subjecti).yhat((si-1)*length(cohs)+ci),...
                        'o','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 750 ms'])
                    hold on
                end
            end
            if includeEarlyInitCohPertTime
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        plot(gainRegression(subjecti).y(3*length(cohs)+(si-1)*length(cohs)+ci),...
                            gainRegression(subjecti).yhat(3*length(cohs)+(si-1)*length(cohs)+ci),...
                            'd','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 50 ms'])
                        hold on
                    end
                end
            end
        elseif dynCohCollate
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi),gainRegression(subjecti).yhat(seqi),...
                    'd','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+length(sequences)),gainRegression(subjecti).yhat(seqi+length(sequences)),...
                    'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+2*length(sequences)),gainRegression(subjecti).yhat(seqi+2*length(sequences)),...
                    's','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                hold on
            end           
        end
        plotUnity;
        axis equal
        axis tight
        xlabel('Behaioral gain')
        ylabel('Estimated gain')
        set(gca,'TickDir','out')
        
        subplot(length(subjects),2,2+2*(subjecti-1))
        stem(regressorClusters,gainRegression(subjecti).B(1:end-1))
        xlabel('Cluster')
        ylabel('Regression weight')
        axis squaregg
        set(gca,'TickDir','out')
    end
    
    %% Topology of gain
    hGainTopology = figure('Name','Gain topology','Position',[1030 669 1312 468]);
    for subjecti = 1:length(subjects)
        subplot(1,length(subjects),subjecti)
        tempGain = nan(size(Y,1),1);
        for typei = 1:NumClusters
            tempGain(idx'==typei) = gainRegression(subjecti).B(typei)/(max(gainRegression(subjecti).B(1:NumClusters)) - min(gainRegression(subjecti).B(1:NumClusters)));
        end
        scatter(Y(:,1),Y(:,2),abs(tempGain)*200,tempGain,'filled')
        hold on
        xlabel('tSNE 1')
        ylabel('tSNE 2')
        axis equal
        axis square
        set(gca,'TickDir','out')
        colormap autumn
        title(subjects{subjecti})
    end
        
        
            
    
    %% Visualize dimensionality reduction
    
    hFunctionalTopology = figure('Name','Functional topography','Position',[1396 220 560 1109]);
    subjectMarkers = {'o','square','diamond','^','<','>','pentagram','hexagram'};
    subplot(4,2,[1,2,3,4])
    for typei = 1:NumClusters
        for subjecti = 1:length(subjects)
            plot(Y(idx'==typei & cellID(:,1,4) == subjecti,1),Y(idx'==typei & cellID(:,1,4) == subjecti,2),...
                subjectMarkers{subjecti},'Color',colorWheel(typei,:),...
                'MarkerFaceColor',colorWheel(typei,:),...
                'DisplayName',[subjects{subjecti} ', cluster ' num2str(typei)]);
            hold on
        end
    end
    axis square
    xlabel('tSNE 1')
    ylabel('tSNE 2')
    
    subplot(4,2,[5,6])
    if initCohCollate
        for i = 1:NumClusters
            plot(initCoh.neuron_t,...
                nanmean(squeeze(Rinit(:,3,3,idx == i))./...
                repmat(max(squeeze(Rinit(:,3,3,idx==i)),[],1),[size(Rinit,1) 1]),2),...
                'Color',colorWheel(i,:),...
                'DisplayName',['InitCoh, speed = 20 deg/s, Coh  = 100%, cluster ' num2str(i)])
            hold on
        end
    end
    title('Response by cluster in initCoh, 20 deg/s target and 100% coherence')
    xlabel('Time from motion onset (ms)')
    ylabel('Normalized response')
    
    subplot(4,2,[7,8])
    if dynCohCollate
        for i = 1:NumClusters
            plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
                nanmean(squeeze(Rdyn(:,5,idx == i))./...
                repmat(max(squeeze(Rdyn(:,5,idx==i)),[],1),[size(Rdyn,1) 1]),2),...
                'Color',colorWheel(i,:),...
                'DisplayName',['DynCoh, speed = 10 deg/s, sequence 5, cluster ' num2str(i)])
            hold on
        end
    end
    title('Response by cluster in dynCoh, control sequence')
    xlabel('Time from motion onset (ms)')
    ylabel('Normalized response')
    
    %% Responses for each speed and coherence in initCoh
    hClusterPSTHsInitCoh = figure('Name','Coherence effect, initCoh','Position',[19 196 3*570 1133]);
    ind = 0;
    for i = 1:NumClusters
        for si = 1:size(Rinit,2)
            ind = ind+1;
            subplot(NumClusters,size(Rinit,2),ind)
            if initCohCollate
                for ci = 1:size(Rinit,3)
                    tempR = nanmean(squeeze(Rinit(:,si,ci,idx == i)),2);
                    plot(initCoh.neuron_t,...
                        tempR,...
                        'Color',initColors(ci,:),...
                        'DisplayName',['Speed = ' num2str(speeds(si)) ' deg/s, Coh = ' num2str(cohs(ci)) '%'])
                    hold on
                end
            end
            axis tight
            lims(ind,:) = axis;
        end
    end
    ind = 0;
    for i = 1:NumClusters
        for si = 1:size(Rinit,2)
            ind = ind+1;
            subplot(NumClusters,size(Rinit,2),ind)
            axis([min(lims(:,1)),max(lims(:,2)),min(lims(:,3)),max(lims(:,4))])
            xlabel('Time from motion onset (ms)')
            ylabel('Normalized response')
        end
    end
    
    %% Responses for each cluster in dynCoh
    hClusterPSTHsDynCoh = figure('Name','Coherence effect, dynCoh','Position',[1956 196 570 1133]);
    for i = 1:NumClusters
        subplot(NumClusters,1,i)
        if dynCohCollate
            for seqi = 1:size(Rdyn,2)
                tempR = nanmean(squeeze(Rdyn(:,seqi,idx == i))./...
                    repmat(max(squeeze(Rdyn(:,seqi,idx==i)),[],1),[size(Rdyn,1) 1]),2);
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
                    tempR,...
                    'Color',colors(seqi,:),'DisplayName',['Sequence ' num2str(seqi)])
                hold on
            end
        end
        axis tight
        lims2(i,:) = axis;
    end
    for i = 1:NumClusters
        subplot(NumClusters,1,i)
        for seqi = 1:size(Rdyn,2)
            axis([min(lims2(:,1)),max(lims2(:,2)),min(lims2(:,3)),max(lims2(:,4))])
            plotVertical([450 750 1050]);
            xlabel('Time from motion onset (ms)')
            ylabel('Normalized response')
        end
    end
    
    %% Validate functional topography
    figure
    
    % select random units
    nExamps = 3;
    unitInds = randsample(size(Rinit,4),nExamps);
    
    % Find K nearest neighbors in functional space
    neighbors = knnsearch(Y,Y(unitInds,:),'K',Kneighbors);
    
    % For each neuron, plot PSTH of nearest neighbors
    for exUnit = 1:nExamps
        subplot(nExamps,2,1+(exUnit-1)*2)
        for ni = 1:Kneighbors
            if dynCohCollate
                zR = (squeeze(Rdyn(:,5,neighbors(exUnit,ni)))-mean(squeeze(Rdyn(:,5,neighbors(exUnit,ni)))))/std(squeeze(Rdyn(:,5,neighbors(exUnit,ni))));
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),zR,'Color',colors(ni,:))
                hold on
            elseif initCohCollate
                zR = (squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))-mean(squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))))/std(squeeze(Rinit(:,3,3,neighbors(exUnit,ni))));
                plot(initCoh.neuron_t(initCoh.neuron_t<=1350),zR(initCoh.neuron_t<=1350),'Color',colors(ni,:))
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel('z-score')
        
        subplot(nExamps,2,2+(exUnit-1)*2)
        plot(Y(:,1),Y(:,2),'k.')
        hold on
        for ni = 1:Kneighbors
            plot(Y(neighbors(exUnit,ni),1),Y(neighbors(exUnit,ni),2),'o','Color',colors(ni,:))
        end
    end
    
    %% topography
    hChamberTopology = figure;
    for subjecti = 1:length(subjects)
        subplot(1,length(subjects),subjecti)
        randScale = 0.08;
        locations2 = locations(cellID(:,1,4) == subjecti,:);
        locations2(locations2(:,1)>1,:) = [locations2(locations2(:,1)>1,2), locations2(locations2(:,1)>1,1), locations2(locations2(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
        locationsRand = locations2 + ...
            [randScale*nanstd(locations2(:,1))*randn(size(locations2,1),1), ...
            randScale*nanstd(locations2(:,2))*randn(size(locations2,1),1), ...
            0*nanstd(locations2(:,3))*randn(size(locations2,1),1)];                  % Add randomness to a-p and m-l locations to make
        for uniti = 1:size(locations2,1)
            plot3(locationsRand(uniti,1),locationsRand(uniti,2),locationsRand(uniti,3)/1000,...
                'o','Color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',colorWheel(idx(uniti),:,:));
            hold on
        end
        grid on
        axis equal
        xlabel('Anterior/postieror (mm)')
        ylabel('Medial/lateral (mm)')
        zlabel('Depth (mm)')
        title([subjects{subjecti} ' recoding locations'])
    end
    
    %% Functional vs physical topography
    figure('Name','Functional vs physical topography','Position',[1153 924 1373 405])
    for subjecti = 1:length(subjects)
        subplot(length(subjects),2,1+2*(subjecti-1))
        samps = randsample(length(euclidLoc{subjecti}),500,false);
        plot(euclidLoc{subjecti}(samps),euclidY{subjecti}(samps),subjectMarkers{subjecti})
        hold on
        ax = axis;
        text((ax(2)-ax(1))*0.8+ax(1),(ax(4)-ax(3))*0.9+ax(3),['R = ' num2str(locCor(1,2))])
        text((ax(2)-ax(1))*0.8+ax(1),(ax(4)-ax(3))*0.85+ax(3),['p = ' num2str(locCorP(1,2))])
        xlabel('Physical distance (mm)')
        ylabel('Functional distance (a.u.)')
        title(['Subject ' subjects{subjecti}])
        
        subplot(length(subjects),2,2+2*(subjecti-1))
        for ni = 1:size(nnIdx{subjecti},1)
            mdist(ni) = mean(pdist(Y(nnIdx{subjecti}(ni,:),:)));
            randIdx = [ni; randsample(size(Y,1),Kneighbors-1)];
            mdistRand(ni) = mean(pdist(Y(randIdx,:)));
        end
        histogram(mdist,linspace(min([mdist,mdistRand]),max([mdist,mdistRand]),50))
        hold on
        histogram(mdistRand,linspace(min([mdist,mdistRand]),max([mdist,mdistRand]),50))
        xlabel('Functional distance')
        ylabel('N')
        legend({[num2str(Kneighbors) ' nearest in chamber'],'Random'})
        title(['Subject ' subjects{subjecti}])
    end
    
    %% Heat map of population firing rates, sorted by embedding location 
    if initCohCollate
        hInitCohRatesHeatMap = figure('Name','Firing rate heat map, initCoh','Position',[864 67 560 1169]);
        cax = [Inf,-Inf];
        RinitNorm = Rinit./max(Rinit,[],[1,2,3]);
        for speedi = 1:length(speeds)
            for cohi = 1:length(cohs)
                subplot(length(speeds),length(cohs),cohi + (speedi-1)*length(cohs))
                Rtemp = squeeze(RinitNorm(:,speedi,cohi,:));
                %             Rtemp = Rtemp./max(Rtemp,[],1);
                imagesc(initCoh.neuron_t,1:size(Rtemp,2),Rtemp(:,thetaSort)')
                xlabel('Time from motion onset (ms)')
                ylabel('Neuron #')
                
                caxTemp = caxis;
                cax(1) = min([caxTemp(1) cax(1)]);
                cax(2) = max([caxTemp(2) cax(2)]);
            end
        end
        for speedi = 1:length(speeds)
            for cohi = 1:length(cohs)
                subplot(length(speeds),length(cohs),cohi + (speedi-1)*length(cohs))
                caxis(cax)
                set(gca,'TickDir','out')
            end
        end
    end
    
    if dynCohCollate
        hDynCohRatesHeatMap = figure('Name','Firing rate heat map, dynCoh','Position',[864 67 560*5/3 1169/3]);
        cax = [Inf,-Inf];
        RdynNorm = Rdyn./max(Rdyn,[],[1,2]);
        for seqi = 1:length(sequences)
            subplot(1,length(sequences),seqi)
            Rtemp = squeeze(RdynNorm(:,seqi,:));
            %             Rtemp = Rtemp./max(Rtemp,[],1);
            imagesc(dynCoh.neuron_t(dynCoh.neuron_t<=1350),1:size(Rtemp,2),Rtemp(:,thetaSort)')
            xlabel('Time from motion onset (ms)')
            ylabel('Neuron #')
            
            caxTemp = caxis;
            cax(1) = min([caxTemp(1) cax(1)]);
            cax(2) = max([caxTemp(2) cax(2)]);
        end
        for seqi = 1:length(sequences)
            subplot(1,length(sequences),seqi)
            caxis(cax)
            set(gca,'TickDir','out')
        end
    end
    
    
    %% Save figures
    if saveFigures
        
        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' [subjects{:}] ...
            '/neuronTyping/' datestr(now,'yyyymmdd')];
        if ~exist(saveLocation,'dir')
            mkdir(saveLocation)
        end
        
        for uniti = 1:size(exUnits,1)
            exIndex = exUnits(uniti,:);
            cellID2 = squeeze(cellID(:,1,1:2));
            listIndex = find(ismember(cellID2, exIndex, 'rows'));
            savefig(hExUnits(uniti),[saveLocation '/exNeuron_' num2str(exIndex(1)) '_' num2str(exIndex(2)) '_' subjects{cellID(listIndex,1,4)} '.fig'])
        end
        for i = 1:NumClusters
            savefig(gvrh(i), [saveLocation '/gainVsClusterFiringRates' num2str(i) '.fig'])
            savefig(hClusterPSTH(i), [saveLocation '/clusterPSTH' num2str(i) '.fig'])
        end
        savefig(hGainDecodingFromClusters , [saveLocation '/gainDecodingClusters.fig'])
        savefig(hFunctionalTopology ,[saveLocation '/functionalTopology.fig'])
        savefig(hGainTopology,[saveLocation '/gainTopology.fig'])
        savefig(hChamberTopology, [saveLocation '/chamberTopology.fig'])
        if initCohCollate
            savefig(hClusterPSTHsInitCoh ,[saveLocation '/clusterPSTHsInitCoh.fig'])
            savefig(hInitCohRatesHeatMap ,[saveLocation '/initCohRatesHeatMap.fig'])
        end
        if dynCohCollate
            savefig(hDynCohRatesHeatMap ,[saveLocation '/dynCohRatesHeatMap.fig'])  
            savefig(hClusterPSTHsDynCoh ,[saveLocation '/clusterPSTHsDynCoh.fig'])
        end
        
    end
end