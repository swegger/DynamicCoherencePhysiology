function neuronTypingAnalysis_initCoh(varargin)
%%
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'binT',0:100:1400)
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'initWin',[150 200])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'win',-200:200)
addParameter(Parser,'shuffleN',100)
addParameter(Parser,'sampleRate',40)
addParameter(Parser,'ClusterMethod','densityClust')
addParameter(Parser,'pertWin',250)
addParameter(Parser,'resultsFile','none')
addParameter(Parser,'testInitGain',true)
addParameter(Parser,'dcpInitCohPertFile','dcpObjectsPertTemp.mat')
addParameter(Parser,'initCohPertFile',[])
addParameter(Parser,'initSpeed',10)
addParameter(Parser,'includeEarlyInitCohPertTime',false)
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'perturbations',[0; 4; 6; 8])
addParameter(Parser,'NumClusters',8)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
binT = Parser.Results.binT;
directions = Parser.Results.directions;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
win = Parser.Results.win;
shuffleN = Parser.Results.shuffleN;
sampleRate = Parser.Results.sampleRate;
ClusterMethod = Parser.Results.ClusterMethod;
pertWin = Parser.Results.pertWin;
resultsFile = Parser.Results.resultsFile;
testInitGain = Parser.Results.testInitGain;
initSpeed = Parser.Results.initSpeed;
includeEarlyInitCohPertTime = Parser.Results.includeEarlyInitCohPertTime;
dcpInitCohPertFile = Parser.Results.dcpInitCohPertFile;
initCohPertFile = Parser.Results.initCohPertFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
perturbations = Parser.Results.perturbations;
NumClusters = Parser.Results.NumClusters;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Load dcp object file
if strcmp(resultsFile,'none')
    load(dcpObjectFile);
    fileExist = false;
else
    fileExist = exist(resultsFile,'file');
    if fileExist
        load(resultsFile)
        fileExist = true;
    else
        error(['File ' resultsFile ' does not exist in path.'])
    end
end

%% Find the gain from behavioral data
if testInitGain
    disp('Determining gain from initiate coherence trials...')
    if ~exist('initGain','var')
        if ~isempty(initCohPertFile) && isfile(initCohPertFile)
            load(initCohPertFile,'init')
        else
            dcpInitCohPert = load(dcpInitCohPertFile);
            [init,~] = initialCohPertBehavioralAnalysis(dcp{1}.sname,'dcp',dcpInitCohPert.dcp(1:end),...
                'outlierReject',false,'win',[150 200],'keep_pert3always0deg',false,'directions',[0,180],'pertWin',pertWin);
        end
        initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
    end
end


%% Results file check start
if ~fileExist
    
    disp('Neural analysis loop...')
    
    %% Get mean and covariance of each unit
    [Rinit, ~, cellID, passCutoff, locations] = collateFiringRates(dcp,...
        'sourceDirectory',sourceDirectory,'directions',directions,'chanMap',chanMap,...
        'rateCutoff',rateCutoff,'checkUnitType',checkUnitType,...
        'initCohCollate',true,'dynCohCollate',false);
    
    temp = load([sourceDirectory '/' dcp{1}.datapath(end-8:end-1) 'obj/initCoh' ...
             dcp{1}.datapath(end-8:end)]);
    initCoh_t = temp.initCoh.neuron_t;
    
%     passCutoff = nan(1000,1);
%     Rinit = nan(1701,3,3,1000);
%     locations = nan(1000,3);
%     cellID = nan(1000,100,3);
%     indx = 1;
%     for filei = 1:length(dcp)
%         disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
%         
%         % Add probe info
%         dcp{filei} = addProbeInfo(dcp{filei});
%         
%         % InitCoh data
%         load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
%             dcp{filei}.datapath(end-8:end)])
%         
%         if ~isnan(rateCutoff)
%             initCoh = findActive(initCoh,rateCutoff,initCoh.cutWindow);
%         end
%         
%         if ~isempty(initCoh.R)
%             
%             passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;
%             
%             
%             e = sqrt(vertcat(initCoh.eye(:).hvel).^2 + vertcat(initCoh.eye(:).vvel).^2)';
%             eInit = nanmean(e(initCoh.eye_t >= initWin(1) & initCoh.eye_t <= initWin(2),:),1);
%             
%             % Get data for each neuron
%             if checkUnitType && isprop(dcp{filei},'unitTypes')
%                 unitInd = find(strcmp(dcp{filei}.unitTypes,'good'));
%             else
%                 unitInd = 1:length(initCoh.preferredDirectionRelative);
%             end
%             for uniti = unitInd
%                 ind = find(dir == initCoh.preferredDirectionRelative(uniti));
%                 Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
%                 
%                 for j = 1:length(initCoh.unitIndex)
%                     cellID(indx,j,1) = filei;
%                     cellID(indx,j,2) = initCoh.unitIndex(uniti);
%                     cellID(indx,j,3) = initCoh.unitIndex(j);
%                 end
%                 
%                 if isempty(initCoh.location)
%                     x = NaN;
%                     y = NaN;
%                     z = NaN;
%                 elseif length(initCoh.location.x)==24
%                     siteIndex = chanMap(initCoh.chansIndex(uniti) == chanMap(:,1),2);
%                     x = initCoh.location.x(siteIndex);
%                     y = initCoh.location.y(siteIndex);
%                     depth = -initCoh.location.depth(siteIndex);
%                 elseif length(initCoh.location.x) > 1
%                     siteIndex = floor(initCoh.chansIndex(uniti)/4)+1;
%                     tempIndex = find(~isnan(initCoh.location.x));
%                     if siteIndex>length(tempIndex)
%                         x = NaN;
%                         y = NaN;
%                         depth = NaN;
%                     else
%                         x = initCoh.location.x(tempIndex(siteIndex));
%                         y = initCoh.location.y(tempIndex(siteIndex));
%                         depth = -initCoh.location.depth(tempIndex(siteIndex));
%                     end
%                 else
%                     x = initCoh.location.x;
%                     y = initCoh.location.y;
%                     depth = -initCoh.location.depth;
%                 end
%                 locations(indx,:) = [x,y,depth];
%                 
%                 indx = indx+1;
%             end
%         end
%     end
%     
%     Rinit = Rinit(:,:,:,1:indx-1);
%     locations = locations(1:indx-1,:);
%     passCutoff = logical(passCutoff(1:indx-1));
%     cellID = cellID(1:indx-1,:,:);
    
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
    
    %% Use nonlinear dimensionality reduction to attempt to identify neuron 'types'
    NumDimensions = 2;
    Perplexity = 10;
    Y = tsne(Rinit3','Distance','cosine','NumDimensions',NumDimensions,'Perplexity',Perplexity);
    
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
    
    %% Distances in functional space vs physical space
    locations2 = locations;
    locations2(locations(:,1)>1,:) = [locations(locations(:,1)>1,2), locations(locations(:,1)>1,1), locations(locations(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
    locations2(:,end) = locations2(:,end)/1000;
    
    euclidLoc = pdist(locations2);
    euclidY = pdist(Y);
    
    % Population level distance correlation
    [locCor,locCorP] = corrcoef(euclidLoc(~isnan(euclidLoc)&~isnan(euclidY)),euclidY(~isnan(euclidLoc)&~isnan(euclidY)));
    
    % K nearest neighbor distances
    Kneighbors = 10;
    nnIdx = knnsearch(locations2,locations2,'K',Kneighbors+1);
    for ni = 1:size(nnIdx,1)
        mdist(ni) = mean(pdist(Y(nnIdx(ni,2:end),:)));
        randIdx = [randsample(size(Y,1),Kneighbors+1)];
        mdistRand(ni) = mean(pdist(Y(randIdx,:)));
    end
    
    %% PCA
    CR = cov(Rinit3);
    [vecs,vals] = eig(CR);
    vecs = fliplr(vecs);
    vals = flipud(diag(vals));
    reconstructionN = 5;
    
    sz = size(Rinit);
    RlowD = nan([sz(1:3),reconstructionN]);
    for si = 1:size(Rinit,2)
        for ci = 1:size(Rinit,3)
            RlowD(:,si,ci,:) = permute(vecs(:,1:reconstructionN)'*squeeze(Rinit2(:,si,ci,:))',[2,3,4,1]);
        end
    end
    
    %% Targeted dimensionality reduction
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(Rinit3);
    zInit = Rinit2./repmat(std(Rinit2,[],[1,2,3]),[size(Rinit2,1),size(Rinit2,2),size(Rinit2,3)]);
    
    % Linear model
    Binit = nan(4,size(zInit,4),size(zInit,1));
    lagT = 100;
    [Cohs Spds] = meshgrid(cohs,speeds);
    if testInitGain
        tempGain = initGain(:,:,3);
    else
        tempGain = zeros(size(Spds));
    end
    for uniti = 1:size(zInit,4)
        for ti = 1:size(zInit,1)
            ztemp = reshape(zInit(ti,:,:,uniti),[numel(Spds),1]);
            Binit(:,uniti,ti) = regress(ztemp,[Spds(:),Cohs(:),tempGain(:),ones(numel(Spds),1)]);
        end        
    end
    
    Dtemp = nan(size(COEFF,1),size(COEFF,1));
    BinitPCA = nan(4,size(Binit,2),size(Binit,3));
    for n = 1:24
        Dtemp(:,:,n) = COEFF(:,n)*COEFF(:,n)';
    end
    D = sum(Dtemp,3);
    for ti = 1:size(zInit,1)
        BinitPCA(:,:,ti) = permute(D*Binit(:,:,ti)',[2,1,3]);
        normB(ti) = norm(BinitPCA(:,:,ti));
    end
    maxNormIndx = find(normB == max(normB));
    BinitMax = BinitPCA(:,:,maxNormIndx);
    BinitTi = BinitPCA(:,:,initCoh_t==750);
    [BinitOrth,~] = qr(BinitTi');
    sz = size(Rinit2);
    Xinit = COEFF(:,1:10)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])';
    pInit = reshape((BinitOrth(:,1:4)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
    
    %% Find targeted dimensions that optimize the representation across time
    
    % Optimized decoder across time
    taccept = initCoh_t >= 150 & initCoh_t <= 1200;
    bOpt = nan(4,size(zInit,4));
    bOptIsolated = nan(2,size(zInit,4));
    for uniti = 1:size(zInit,4)
        temp = reshape(permute(zInit(taccept,:,:,uniti),[2,3,1]),[numel(Spds)*sum(taccept),1]);
        bOpt(:,uniti) = regress(temp,...
            repmat([Spds(:),Cohs(:),tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
        bOptIsolated(:,uniti) = regress(temp,...
            repmat([tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
    end
        
    % Find unit vector 
    normBopt = nan(size(bOpt));
    normBoptIsolated = nan(size(bOptIsolated));
    for dimi = 1:size(bOpt,1)
        normBopt(dimi,:) = bOpt(dimi,:)/norm(bOpt(dimi,:));
    end
    normBoptIsolated(1,:) = bOptIsolated(1,:)/norm(bOptIsolated(1,:));
    normBoptIsolated(2,:) = bOptIsolated(2,:)/norm(bOptIsolated(2,:));
    
    % Project into major PC dimensions
    bOptPCA = permute(D*bOpt',[2,1]);
    
    % Find orthonormal basis
    [bOptOrth,~] = qr(bOpt');
    
    % Project along optimal subspace dimensions
    sz = size(zInit);
    pBopt = reshape((bOptOrth(:,1:4)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
    pBoptIsolated = reshape((normBoptIsolated*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),2]);
    
    
end

%% Save results file
if saveResults
    
    saveLocation = ['/home/seth/Projects/DynamicCoherencePhysiology/' dcp{1}.sname ...
        '/initCohDynCohComp'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/neuronTyping_initCoh' datestr(now,'yyyymmdd')],'-v7.3')
    
end

%% Simple avaraging
a = Rinit(initCoh_t<=1350,:,:,:);

figure('Name','Grand average','Position',[486 733 516 420])
minMax = [Inf,0];
for speedi = 1:length(speeds)
    subplot(length(speeds),1,speedi)
    for cohi = 1:length(cohs)
        plot(initCoh_t(initCoh_t<=1350),nanmean(squeeze(a(:,speedi,cohi,:)),2)*1000,...
            'Color',[speedColors(speedi,:), cohs(cohi)/100])
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

speedi = find(speeds == 10);

gvmrh = figure('Name',['Grand average']);
for si = 1:length(speeds)
    initRatesTemp = squeeze(nanmean(a(initCoh_t >= 700 & initCoh_t <= 800,si,:,:),[1,4]))*1000;
    plot(squeeze(initGain(si,:,3)),initRatesTemp,...
        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    hold on
end
if includeEarlyInitCohPertTime
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(initCoh_t >= 100 & initCoh_t <= 200,si,:,:),[1,4]))*1000;
        plot(squeeze(initGain(si,:,2)),initRatesTemp,...
            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    end
end

xlabel('Behavioral gain (unitless)')
ylabel('Spikes/s')


%% Simple avaraging, by cluster
gainRegression.y = [];
gainRegression.x = nan(length(speeds)*length(cohs),NumClusters);
for i = 1:NumClusters
    a = Rinit(initCoh_t<=1350,:,:,idx == i);
    
    figure('Name',['Cluster ' num2str(i)])    
    subplot(2,1,1)
    plot(initCoh_t(initCoh_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)),2)*1000,'Color',initColors(1,:))
    hold on
    plot(initCoh_t(initCoh_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
    plot(initCoh_t(initCoh_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)),2)*1000,'Color',initColors(3,:))
    xlabel('Time from motion onset (ms)')
    ylabel('Spike/s')

    subplot(2,1,2)
    for si = 1:length(speeds)
        plot(initCoh_t(initCoh_t<=1350),nanmean(squeeze(a(:,si,2,:)),2)*1000,'Color',speedColors(si,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    
    gvrh = figure('Name',['Cluster ' num2str(i)]);
    for si = 1:length(speeds)
        initRatesTemp = nanmean(squeeze(nanmean(a(initCoh_t >= 700 & initCoh_t <= 800,si,:,:),1)),2)*1000;
        gainRegression.x((si-1)*length(cohs)+1:si*length(cohs),i) = initRatesTemp;
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        hold on
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            initRatesTemp = nanmean(squeeze(nanmean(a(initCoh_t >= 100 & initCoh_t <= 200,si,:,:),1)),2)*1000;
            gainRegression.x(3*length(cohs) + (si-1)*length(cohs)+1:3*length(cohs)+si*length(cohs),i) = initRatesTemp;
            plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
    end
    xlabel('Behavioral gain (unitless)')
    ylabel('Spikes/s')
end

%% Gain decoding performance
regressorClusters = 1:NumClusters; %[2,7,8];
gainRegression.x_mean = mean(gainRegression.x,1);
gainRegression.x_std = std(gainRegression.x,[],1);
gainRegression.z = gainRegression.x;
if testInitGain
    for si = 1:length(speeds)
        gainRegression.y(1+(si-1)*length(cohs):si*length(cohs)) = initGain(si,:,3);
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            gainRegression.y(3*length(cohs)+(si-1)*length(cohs)+1:3*length(cohs)+si*length(cohs)) = initGain(si,:,2);
        end
    end
end
gainRegression.B = regress(gainRegression.y(:),[gainRegression.z(:,regressorClusters) ones(size(gainRegression.x,1),1)]);
gainRegression.yhat = [gainRegression.z(:,regressorClusters) ones(size(gainRegression.x,1),1)]*gainRegression.B;
gainRegression.Bnonneg = lsqnonneg(gainRegression.z,gainRegression.y');

figure('Name','Gain decoding','Position',[1633 927 888 395])
subplot(1,2,1)
if testInitGain
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            plot(gainRegression.y((si-1)*length(cohs)+ci),...
                gainRegression.yhat((si-1)*length(cohs)+ci),...
                'o','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
            hold on
        end
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            for ci = 1:length(cohs)
                plot(gainRegression.y(3*length(cohs)+(si-1)*length(cohs)+ci),...
                    gainRegression.yhat(3*length(cohs)+(si-1)*length(cohs)+ci),...
                    'd','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
                hold on
            end
        end
    end
end
plotUnity;
axis equal
axis tight
xlabel('Behaioral gain')
ylabel('Estimated gain')

subplot(1,2,2)
stem(regressorClusters,gainRegression.B(1:end-1))
xlabel('Cluster')
ylabel('Regression weight')
axis square

%% PCA plotting

figure
for ri = 1:reconstructionN
    for si = 1:size(RlowD,2)
        subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
        for ci = 1:size(RlowD,3)
            plot(initCoh_t,RlowD(:,si,ci,ri),'Color',initColors(ci,:))
            hold on
        end
        templims(si,:) = axis;
    end

    for si = 1:size(RlowD,2)
        subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
        ylim([min(templims(:,3)) max(templims(:,4))]);
    end
        
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
end

figure
for si = 1:size(RlowD,2)
    for ci = 1:size(RlowD,3)
        plot3(RlowD(:,si,ci,1),RlowD(:,si,ci,2),RlowD(:,si,ci,3),'-','Color',[speedColors(si,:) cohs(ci)/100])
        hold on
        scatter3(RlowD(1,si,ci,1),RlowD(1,si,ci,2),RlowD(1,si,ci,3),100,speedColors(si,:),...
            'filled','o','MarkerEdgeAlpha',cohs(ci)/100,'MarkerFaceAlpha',cohs(ci)/100)
        scatter3(RlowD(end,si,ci,1),RlowD(end,si,ci,2),RlowD(end,si,ci,3),100,speedColors(si,:),...
            'filled','s','MarkerEdgeAlpha',cohs(ci)/100,'MarkerFaceAlpha',cohs(ci)/100)
    end
end
grid on
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
view([-75 60])

%% Plot targeted dimensionality reduction
figure('Name','Input related activity','Position',[63 169 1606 1079])
subplot(3,1,1)
for speedi = 1:3
    for cohi = 1:3
        plot(initCoh_t,pInit(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
axis tight
ax(1,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Speed related activity (a.u.)')

subplot(3,1,2)
for cohi = 1:3
    for speedi = 1:length(speeds)
        plot(initCoh_t,pInit(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
axis tight
ax(1,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Coherence related activity (a.u.)')

subplot(3,1,3)
for cohi = 1:3
    for speedi = 1:length(speeds)
        plot(initCoh_t,pInit(:,speedi,cohi,3),'-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
axis tight
ax(1,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Gain related activity (a.u.)')

%% Targeted dimension activity vs gain
dimNames = {'Speed','Coherence','Gain','Offset'};
gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
for targetedDim = 1:3
    tempData = [];
    subplot(3,1,targetedDim)
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(pInit(initCoh_t >= 700 & initCoh_t <= 800,si,:,targetedDim),1));
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        hold on
%         if si == 3
%             tempData = [tempData; initGain(si,1:2,3)',initRatesTemp(1:2)];
%         else
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%         end
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(initCoh_t >= 100 & initCoh_t <= 200,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
        end
    end
    gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
    axis square
    ax = axis;
    text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
    xlabel('Behavioral gain (unitless)')
    ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
end

%% Plot activity along orthonormal subspace defined by targeted optimized across time points
figure('Name','Activity along targeted dimensions found by optimizing over time',...
    'Position',[73 164 1616 1072])
dimensionLabels = {'Speed','Coherence','Gain'};
for di = 1:3
    optH(di) = subplot(3,1,di);
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            plot(initCoh_t,pBopt(:,si,ci,di),'-','Color',[speedColors(si,:) cohs(ci)/100])
            hold on
        end
    end
    axis tight
    xlabel('Time from motion onset (ms)')
    ylabel('Activity along targeted dimension')
    title(dimensionLabels{di})
    
end

linkaxes(optH,'xy')

%% Plot behavioral gain vs activity along targeted dimensions
dimNames = {'Speed','Coherence','Gain','Offset'};
gvthOpt = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
for targetedDim = 1:3
    tempData = [];
    subplot(3,1,targetedDim)
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(pBopt(initCoh_t >= 700 & initCoh_t <= 800,si,:,targetedDim),1));
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        hold on
        tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pBopt(initCoh_t >= 100 & initCoh_t <= 200,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
        end
    end
        
    gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
    axis square
    ax = axis;
    text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
    xlabel('Behavioral gain (unitless)')
    ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
end

%% Visualize dimensionality reduction

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

vecInds = [1,2,3];
% vColor = sqrt(vecs(:,vecInds).^2./repmat(max(vecs(:,vecInds).^2,[],2),[1,3]));

figure('Name','Functional topography','Position',[1396 220 560 1109]);
if NumDimensions > 2
    dimsTemp = [1,2,3];
    subplot(2,3,[1,2,4,5])
%     scatter3(Y(:,dimsTemp(1)),Y(:,dimsTemp(2)),Y(:,dimsTemp(3)),50,vColor,'filled')
    hold on
%     scatter3(Y(:,dimsTemp(1)),Y(:,dimsTemp(2)),Y(:,dimsTemp(3)),50,colors(idx,:))
    axis square
    xlabel('tSNE 1')
    ylabel('tSNE 2')

    subplot(2,3,[3,6])
    for i = 1:NumClusters
        plot(initCoh_t,...
            nanmean(squeeze(Rinit(:,3,3,idx == i))./...
            repmat(max(squeeze(Rinit(:,3,3,idx==i)),[],1),[size(Rinit,1) 1]),2))
        hold on
    end

else
    subplot(4,2,[1,2,3,4])
    for typei = 1:NumClusters
        plot(Y(idx==typei,1),Y(idx==typei,2),'o','Color',colorWheel(typei,:),...
            'MarkerFaceColor',colorWheel(typei,:));
        hold on
    end
%     scatter(Y(:,1),Y(:,2),50,colors(idx,:),'filled')
%     hold on
%     scatter(Y(:,1),Y(:,2),25,vColor,'filled')
    axis square
    xlabel('tSNE 1')
    ylabel('tSNE 2')

    subplot(4,2,[5,6])
    for i = 1:NumClusters
        plot(initCoh_t,...
            nanmean(squeeze(Rinit(:,3,3,idx == i))./...
            repmat(max(squeeze(Rinit(:,3,3,idx==i)),[],1),[size(Rinit,1) 1]),2),...
            'Color',colorWheel(i,:))
        hold on
    end

end

figure('Name','Coherence effect, initCoh','Position',[19 196 3*570 1133])
ind = 0;
for i = 1:NumClusters
    for si = 1:size(Rinit,2)
        ind = ind+1;
        subplot(NumClusters,size(Rinit,2),ind)
        for ci = 1:size(Rinit,3)
            tempR = nanmean(squeeze(Rinit(:,si,ci,idx == i))./...
                repmat(max(squeeze(Rinit(:,si,ci,idx==i)),[],1),[size(Rinit,1) 1]),2);
            plot(initCoh_t,...
                tempR,...
                'Color',initColors(ci,:))
            hold on
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
        ylabel('z-score')
    end
end

figure('Name','Mean Cluster PSTHs','Position',[1956 219 570 1110])
for i = 1:NumClusters
    subplot(NumClusters,1,i)
%     plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
%         nanmean(squeeze(Rdyn(:,5,idx == i))./...
%         repmat(max(squeeze(Rdyn(:,5,idx==i)),[],1),[size(Rdyn,1) 1]),2),...
%         'Color',colorWheel(i,:))
    plot(initCoh_t(initCoh_t<=1350),...
        nanmean((squeeze(Rinit(initCoh_t<=1350,3,3,idx == i))-nanmean(squeeze(Rinit(initCoh_t<=1350,2,3,idx == i)),1))./...
        repmat(nanstd(squeeze(Rinit(initCoh_t<=1350,2,3,idx==i)),[],1),[size(Rinit(initCoh_t<=1350,:,:,:),1) 1]),2),...
        'Color',colorWheel(i,:))
    xlabel('Time from motion onset (ms)')
    ylabel('z-score')
end

%% Plot zscores by 1D projection of functional topology
thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
[~,thetaSort] = sort(thetas);

figure
imagesc(initCoh_t,1:size(zInit,3),squeeze(zInit(:,2,2,thetaSort))')
xlabel('Time from motion onset (ms)')
ylabel('Sorted neuron #')


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
        zR = (squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))-mean(squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))))/std(squeeze(Rinit(:,3,3,neighbors(exUnit,ni))));
        plot(initCoh_t,zR,'Color',speedColors(3,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel('z-score')

    subplot(nExamps,2,2+(exUnit-1)*2)
    plot(Y(:,1),Y(:,2),'k.')
    hold on
    for ni = 1:Kneighbors
        plot(Y(neighbors(exUnit,ni),1),Y(neighbors(exUnit,ni),2),'o','Color',speedColors(3,:))
    end
end

%% topography
figure;
randScale = 0.08;
vecInds = [1,2,3];
locations2 = locations;
locations2(locations(:,1)>1,:) = [locations(locations(:,1)>1,2), locations(locations(:,1)>1,1), locations(locations(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
locationsRand = locations2 + ...
    [randScale*nanstd(locations2(:,1))*randn(size(locations,1),1), ...
    randScale*nanstd(locations2(:,2))*randn(size(locations,1),1), ...
    0*nanstd(locations2(:,3))*randn(size(locations,1),1)];                  % Add randomness to a-p and m-l locations to make
for uniti = 1:size(locations,1)
    plot3(locationsRand(uniti,1),locationsRand(uniti,2),locationsRand(uniti,3)/1000,...
        'o','Color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',colorWheel(idx(uniti),:,:));
    hold on
end
grid on
axis equal
xlabel('Anterior/postieror (mm)')
ylabel('Medial/lateral (mm)')
zlabel('Depth (mm)')

%% Functional vs physical topography
figure('Name','Functional vs physical topography','Position',[1153 924 1373 405])
subplot(1,2,1)
samps = randsample(length(euclidLoc(:)),500,false);
plot(euclidLoc(samps),euclidY(samps),'ko')
ax = axis;
text((ax(2)-ax(1))*0.8+ax(1),(ax(4)-ax(3))*0.9+ax(3),['R = ' num2str(locCor(1,2))])
text((ax(2)-ax(1))*0.8+ax(1),(ax(4)-ax(3))*0.85+ax(3),['p = ' num2str(locCorP(1,2))])
xlabel('Physical distance (mm)')
ylabel('Functional distance (a.u.)')

subplot(1,2,2)
for ni = 1:size(nnIdx,1)
    mdist(ni) = mean(pdist(Y(nnIdx(ni,:),:)));
    randIdx = [ni; randsample(size(Y,1),Kneighbors-1)];
    mdistRand(ni) = mean(pdist(Y(randIdx,:)));
end
histogram(mdist,linspace(min([mdist,mdistRand]),max([mdist,mdistRand]),50))
hold on
histogram(mdistRand,linspace(min([mdist,mdistRand]),max([mdist,mdistRand]),50))
xlabel('Functional distance')
ylabel('N')
legend({[num2str(Kneighbors) ' nearest in chamber'],'Random'})
