function neuronTypingAnalysis(varargin)
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
addParameter(Parser,'dir',[0 180])
addParameter(Parser,'initWin',[150 200])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'win',-200:200)
addParameter(Parser,'calcCC',false)
addParameter(Parser,'calcRSC',true)
addParameter(Parser,'shuffleN',100)
addParameter(Parser,'sampleRate',40)
addParameter(Parser,'ClusterMethod','densityClust')
addParameter(Parser,'pertWin',250)
addParameter(Parser,'resultsFile','none')
addParameter(Parser,'testInitGain',true)
addParameter(Parser,'dcpInitCohPertFile','dcpObjectsPertTemp.mat')
addParameter(Parser,'initCohPertFile',[])
addParameter(Parser,'dcpDynCohFile',[])
addParameter(Parser,'initSpeed',10)
addParameter(Parser,'includeEarlyInitCohPertTime',false)
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'sequences',[1; 2; 3; 4; 5])
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
dir = Parser.Results.dir;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
win = Parser.Results.win;
calcCC = Parser.Results.calcCC;
calcRSC = Parser.Results.calcRSC;
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
dcpDynCohFile = Parser.Results.dcpDynCohFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
sequences = Parser.Results.sequences;
perturbations = Parser.Results.perturbations;
NumClusters = Parser.Results.NumClusters;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
figure;
colors = colormap('lines');
close(gcf)
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

if ~exist('gain','var') || ~exist('dyn','var')
    disp('Determining gain from dynamic coherence trials...')
    if isempty(dcpDynCohFile)
        [dyn,gain] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcp,'directions',[0,180],'pertWin',pertWin);
    else
        dcpDynCoh = load(dcpDynCohFile);
        [dyn,gain] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcpDynCoh.dcp,'directions',[0,180],'pertWin',pertWin);
    end
    eye_tDyn = dyn.t;
    meanEyeSpeedDyn = dyn.eye.mean;
end

if testInitGain
    disp('Determining gain from initiate coherence trials...')
    if ~exist('initGain','var') || ~exist('init','var')
        if ~isempty(initCohPertFile) && isfile(initCohPertFile)
            load(initCohPertFile,'init')
        else
            dcpInitCohPert = load(dcpInitCohPertFile);
            [init,~] = initialCohPertBehavioralAnalysis(dcp{1}.sname,'dcp',dcpInitCohPert.dcp(1:end),...
                'outlierReject',false,'win',[150 200],'keep_pert3always0deg',false,'directions',[0,180],'pertWin',pertWin);
        end
        eye_t = init.t;
        meanEyeSpeed = init.eye.mean;
        initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
    end
end

%% Results file check start
if ~fileExist
    
    disp('Neural analysis loop...')
    
    %% Get mean and covariance of each unit
    passCutoff = nan(1000,1);
    Rinit = nan(1701,3,3,1000);
    Rdyn = nan(1701,5,1000);
    locations = nan(1000,3);
    cellID = nan(1000,100,3);
    RSCinit = nan(1000,100);
    RSCinit_pval = nan(1000,100);
    RSCdyn = nan(1000,100);
    RSCdyn_pval = nan(1000,100);
    indx = 1;
    for filei = 1:length(dcp)
        disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
        
        % Add probe info
        dcp{filei} = addProbeInfo(dcp{filei});
        
        % InitCoh data
        load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
            dcp{filei}.datapath(end-8:end)])
        
        % DynCoh data
        load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{filei}.datapath(end-8:end)])
        if ~isnan(rateCutoff)
            initCoh = findActive(initCoh,rateCutoff,initCoh.cutWindow);
            dynCoh = findActive(dynCoh,rateCutoff,dynCoh.cutWindow);
        end
        if ~isempty(initCoh.R) && (~isempty(dynCoh.R) || ~any(isnan(dynCoh.R(:))) && size(dynCoh.R,2) > 1)
            if length(initCoh.unitIndex) == length(dynCoh.unitIndex)
                
                passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff | dynCoh.passCutoff;
                
                
                e = sqrt(vertcat(initCoh.eye(:).hvel).^2 + vertcat(initCoh.eye(:).vvel).^2)';
                eInit = nanmean(e(initCoh.eye_t >= initWin(1) & initCoh.eye_t <= initWin(2),:),1);
                
                
                
                if calcRSC
                    [rsc,pval] = spikeCountCorrelationWin(initCoh,'win',[150,450]);
                    [rsc2,pval2] = spikeCountCorrelationWin(dynCoh,'win',[150,450]);
                end
                
                if calcCC
                    [cc, shufflecc] = correlograms(initCoh,sampleRate,initCoh.unitIndex,win,shuffleN);
                    counts = spikeCount(initCoh,initCoh.unitIndex);
                end
                
                % Get data for each neuron
                if checkUnitType && isprop(dcp{filei},'unitTypes')
                    unitInd = find(strcmp(dcp{filei}.unitTypes,'good'));
                else
                    unitInd = 1:length(initCoh.preferredDirectionRelative);
                end
                for uniti = unitInd
                    ind = find(dir == initCoh.preferredDirectionRelative(uniti));
                    Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
                    
                    ind = find(dir == dynCoh.preferredDirectionRelative(uniti));
                    Rdyn(:,:,indx) = dynCoh.R(:,:,uniti,ind);
                    
                    
                    for j = 1:length(initCoh.unitIndex)
                        cellID(indx,j,1) = filei;
                        cellID(indx,j,2) = initCoh.unitIndex(uniti);
                        cellID(indx,j,3) = initCoh.unitIndex(j);
                    end
                    
                    if calcRSC
                        for j = 1:length(initCoh.unitIndex)
                            RSCinit(indx,j) = rsc(uniti,j);
                            RSCinit_pval(indx,j) = pval(uniti,j);
                            RSCdyn(indx,j) = rsc2(uniti,j);
                            RSCdym_pval(indx,j) = pval2(uniti,j);
                        end
                    end
                    
                    if calcCC
                        for j = 1:length(initCoh.unitIndex)
                            CC(indx,j,:) = (cc(:,uniti,j)-mean(shufflecc(:,uniti,j,:),4))/counts(uniti);                            % Shuffle corrected correlogram
                            shuffletemp = reshape(shufflecc(:,uniti,j,:),[numel(shufflecc(:,uniti,j,:)),1]);
                            CC_CI(indx,j,:) = [-1.96*std(shuffletemp),1.96*std(shuffletemp)]/counts(uniti);   % 95% confidence intervals of correlogram
                        end
                    end
                    
                    if isempty(initCoh.location)
                        x = NaN;
                        y = NaN;
                        z = NaN;
                    elseif length(initCoh.location.x)==24
                        siteIndex = chanMap(initCoh.chansIndex(uniti) == chanMap(:,1),2);
                        x = initCoh.location.x(siteIndex);
                        y = initCoh.location.y(siteIndex);
                        depth = -initCoh.location.depth(siteIndex);
                    elseif length(initCoh.location.x) > 1
                        siteIndex = floor(initCoh.chansIndex(uniti)/4)+1;
                        tempIndex = find(~isnan(initCoh.location.x));
                        if siteIndex>length(tempIndex)
                            x = NaN;
                            y = NaN;
                            depth = NaN;
                        else
                            x = initCoh.location.x(tempIndex(siteIndex));
                            y = initCoh.location.y(tempIndex(siteIndex));
                            depth = -initCoh.location.depth(tempIndex(siteIndex));
                        end
                    else
                        x = initCoh.location.x;
                        y = initCoh.location.y;
                        depth = -initCoh.location.depth;
                    end
                    locations(indx,:) = [x,y,depth];
                    
                    indx = indx+1;
                end
            end
        end
    end
    
    Rinit = Rinit(:,:,:,1:indx-1);
    Rdyn = Rdyn(:,:,1:indx-1);
    locations = locations(1:indx-1,:);
    passCutoff = logical(passCutoff(1:indx-1));
    cellID = cellID(1:indx-1,:,:);
    RSCinit = RSCinit(1:indx-1,:);
    RSCinit_pval = RSCinit_pval(1:indx-1,:);
    RSCdyn = RSCdyn(1:indx-1,:);
    RSCdyn_pval = RSCdyn_pval(1:indx-1,:);
    
    %% Remove data that doesn't pass cutoff
    Rinit = Rinit(:,:,:,passCutoff);
    Rdyn = Rdyn(:,:,passCutoff);
    locations = locations(passCutoff,:);
    cellID = cellID(passCutoff,:,:);
    RSCinit = RSCinit(logical(passCutoff),:);
    RSCinit_pval = RSCinit_pval(logical(passCutoff),:);
    RSCdyn = RSCdyn(logical(passCutoff),:);
    RSCdyn_pval = RSCdyn_pval(logical(passCutoff),:);
    if calcCC
        CC = CC(logical(passCutoff),:,:);
        CC_CI = CC_CI(logical(passCutoff),:,:);
    end
    
    %% Remove outlier rates
    m = squeeze(max(Rinit,[],[1,2,3]))*1000;
    m2 = squeeze(max(Rdyn,[],[1,2]))*1000;
    Rinit = Rinit(:,:,:,m<=150 & m2<=150);
    Rdyn = Rdyn(:,:,m<=150 & m2<=150);
    locations = locations(m<=150 & m2<=150,:);
    cellID = cellID(m<=150 & m2<150,:,:);
    RSCinit = RSCinit(m<=150 & m2<150,:);
    RSCinit_pval = RSCinit_pval(m<=150 & m2<=150,:);
    RSCdyn = RSCdyn(m<=150 & m2<=150,:);
    RSCdyn_pval = RSCdyn_pval(m<=150 & m2<=150,:);
    if calcCC
        CC = CC(m<=150 & m2<=150,:,:);
        CC_CI = CC_CI(m<=150 & m2<=150,:,:);
    end
    
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
    
    
    %% Use nonlinear dimensionality reduction to attempt to identify neuron 'types'
    NumDimensions = 2;
    Perplexity = 10;
    Y = tsne([Rinit3; Rdyn3]','Distance','cosine','NumDimensions',NumDimensions,'Perplexity',Perplexity);
    
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
    
    %% calculate Rsc by cluster
    if calcRSC
        RSC2 = RSCinit;
        %     RSC2(RSCinit_pval>0.05) = NaN;
        RSC2(RSC2 == 1) = NaN;
        [TEMP1 TEMP2] = meshgrid(1:NumClusters);
        typesAll = [TEMP1(:) TEMP2(:)];
        RSC_typeConditioned = cell(1,size(typesAll,1));
        mRSC = nan(1,size(typesAll,1));
        for compi = 1:size(typesAll,1)
            types = typesAll(compi,:);%[3,3];
            type1 = find(idx == types(1));
            
            tempIDs = nan(length(type1),1);
            ind = 0;
            for uniti = 1:length(type1)
                filei = cellID(type1(uniti),1,1);
                unitIDs = cellID(type1(uniti),:,3);
                unitIDs = unitIDs(~isnan(unitIDs));
                for unitj = 1:length(unitIDs)
                    unit2 = find(cellID(:,1,1) == filei & cellID(:,1,2) ==  unitIDs(unitj));
                    type2 = idx(unit2);
                    if type2 == types(2)
                        ind = ind+1;
                        tempIDs(ind) = RSC2(type1(uniti),unitj);
                    end
                end
            end
            tempIDs = tempIDs(1:ind);
            RSC_typeConditioned{compi} = tempIDs;
            mRSC(compi) = nanmean(tempIDs);
        end
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
    
    sz = size(Rdyn);
    RDynlowD = nan([sz(1:2),reconstructionN]);
    for seqi = 1:size(Rdyn,2)
        
        RDynlowD(:,seqi,:) = permute(vecs(:,1:reconstructionN)'*squeeze(Rdyn2(:,seqi,:))',[2,3,4,1]);
        
    end
    
    %% Targeted dimensionality reduction
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(Rinit3);
    zInit = Rinit2./repmat(std(Rinit2,[],[1,2,3]),[size(Rinit2,1),size(Rinit2,2),size(Rinit2,3)]);
    zDyn = Rdyn2./repmat(std(Rdyn2,[],[1,2]),[size(Rdyn2,1),size(Rdyn2,2)]);
    dynCohT = dynCoh.coh;
    dynCohT = [60*ones(sum(dynCoh.neuron_t<0),5);dynCohT];
    % Linear model
    Binit = nan(4,size(zInit,4),size(zInit,1));
    Bdyn = nan(3,size(zInit,4),size(zDyn,1));
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
        
        for ti = (lagT+1):size(zDyn,1)
            ztemp = reshape(zDyn(ti,:,uniti),[size(zDyn,2),1]);
            if length(unique(dynCohT(ti-lagT,:))) == 1
                Bdyn(1,uniti,ti-lagT) = NaN;
                Bdyn(2:3,uniti,ti-lagT) = regress(ztemp,[gain(:,3),ones(numel(ztemp),1)]);
            else
                warning('off','all')
                Bdyn(:,uniti,ti-lagT) = regress(ztemp,[dynCohT(ti-lagT,:)',gain(:,3),ones(numel(ztemp),1)]);
                warning('on','all') 
            end
        end
        
    end
    
    Dtemp = nan(size(COEFF,1),size(COEFF,1));
    BinitPCA = nan(4,size(Binit,2),size(Binit,3));
    BdynPCA = nan(3,size(Bdyn,2),size(Bdyn,3));
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
    BinitTi = BinitPCA(:,:,initCoh.neuron_t==750);
    [BinitOrth,~] = qr(BinitTi');
    sz = size(Rinit2);
    Xinit = COEFF(:,1:10)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])';
    pInit = reshape((BinitOrth(:,1:4)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
    
    for ti = 1:size(zDyn,1)
        BdynPCA(:,:,ti) = permute(D*Bdyn(:,:,ti)',[2,1,3]);
        %     normB(ti) = norm(BdynPCA(:,:,ti));
    end
    BdynTi = BdynPCA(:,:,dynCoh.neuron_t==750);
    [BdynOrth,~] = qr(BdynTi');
    sz = size(Rdyn2);
    Xdyn = COEFF(:,1:10)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])';
    pDyn = reshape((BdynOrth(:,1:3)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),3]);
    pDynCross = reshape((BinitOrth(:,1:4)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),4]);
    
    %% Regress targeted projection with coherence
    for ti = 1:size(pInit,1)
        for di = 1:2
            betaInit(:,ti,di) = regress(permute(pInit(ti,2,:,di),[3,1,2,4]),[cohs ones(size(cohs))]);
        end
    end
    betaDyn = nan(2,size(pDynCross,1),2);
    startTi = find(dynCoh.neuron_t == 550);
    endTi = find(dynCoh.neuron_t == 1150);
    for ti = startTi:endTi
        for di = 1:2
            betaDyn(:,ti,di) = regress(permute(pDynCross(ti,:,di),[2,1,3]),[dynCohT(ti-100,:)' ones(size(dynCoh.coh(ti-100,:)'))]);
        end
    end
    
    %% Find targeted dimensions that optimize the representation across time
    
    % Optimized decoder across time
    taccept = initCoh.neuron_t >= 150 & initCoh.neuron_t <= 1200;
    bOpt = nan(4,size(zInit,4));
    bOptIsolated = nan(2,size(zInit,4));
    bOptDynCoh = nan(3,size(zDyn,3));
    for uniti = 1:size(zInit,4)
        temp = reshape(permute(zInit(taccept,:,:,uniti),[2,3,1]),[numel(Spds)*sum(taccept),1]);
        bOpt(:,uniti) = regress(temp,...
            repmat([Spds(:),Cohs(:),tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
        bOptIsolated(:,uniti) = regress(temp,...
            repmat([tempGain(:),ones(numel(Spds),1)],[sum(taccept),1]));
    end
    
    taccept = taccept(dynCoh.neuron_t <=1350);
    for uniti = 1:size(zDyn,3)
        temp = reshape(permute(zDyn(taccept,:,uniti),[2,1]),[size(zDyn,2)*sum(taccept),1]);
        bOptDynCoh(:,uniti) = regress(temp,...
            [reshape(permute(dynCoh.coh(taccept,:),[2,1]),[size(zDyn,2)*sum(taccept),1]) repmat(gain(:,3),[sum(taccept),1]) ones(sum(taccept)*size(zDyn,2),1)]);
    end
    
    
    % Find unit vector 
    normBopt = nan(size(bOpt));
    normBoptDynCoh = nan(size(bOptDynCoh));
    normBoptIsolated = nan(size(bOptIsolated));
    for dimi = 1:size(bOpt,1)
        normBopt(dimi,:) = bOpt(dimi,:)/norm(bOpt(dimi,:));
    end
    for dimi = 1:size(bOptDynCoh,1)
        normBoptDynCoh(dimi,:) = bOptDynCoh(dimi,:)/norm(bOptDynCoh(dimi,:));
    end
    normBoptIsolated(1,:) = bOptIsolated(1,:)/norm(bOptIsolated(1,:));
    normBoptIsolated(2,:) = bOptIsolated(2,:)/norm(bOptIsolated(2,:));
    
    % Project into major PC dimensions
    bOptPCA = permute(D*bOpt',[2,1]);
    bOptDynCohPCA = permute(D*bOptDynCoh',[2,1]);
    
    % Find orthonormal basis
    [bOptOrth,~] = qr(bOpt');
    [bOptDynCohOrth,~] = qr(bOptDynCoh');
    
    % Project along optimal subspace dimensions
    sz = size(zInit);
    pBopt = reshape((bOptOrth(:,1:4)'*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
    pBoptIsolated = reshape((normBoptIsolated*reshape(zInit,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),2]);
    
    sz = size(zDyn);
    pOptDynCross = reshape((bOptOrth(:,1:4)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),4]);
    pOptDyn = reshape((bOptDynCohOrth(:,1:3)'*reshape(zDyn,[prod(sz(1:2)),prod(sz(end))])')',[sz(1),sz(2),3]);
    
    %% Cross correllogram analysis
    if calcCC
        nSig = nansum(CC(:,:,151:200) > repmat(CC_CI(:,:,2),[1,1,50]) | CC(:,:,151:200) < repmat(CC_CI(:,:,1),[1,1,50]),3) + ...
            nansum(CC(:,:,202:251) > repmat(CC_CI(:,:,2),[1,1,50]) | CC(:,:,202:251) < repmat(CC_CI(:,:,1),[1,1,50]),3);
        
        nSig(cellID(:,:,2)==cellID(:,:,3)) = NaN;
        
        nSig2 = nan(size(nSig,1),size(nSig,1));
        connectivity = nan(size(nSig,1),size(nSig,2),4);
        chanN = nan(size(nSig,1),size(nSig,1));
        fileID = nan(size(nSig,1),size(nSig,1));
        for i = 1:size(nSig,1)
            for j = 1:size(nSig,2)
                newi = i;
                newj = find(cellID(:,1,1) == cellID(i,j,1) & cellID(:,1,2) == cellID(i,j,3));
                nSig2(newi,newj) = nSig(i,j);
                
                connectivity(newi,newj,1) = nansum(CC(i,j,201-20:200) > repmat(CC_CI(i,j,2),[1,1,length(201-20:200)]));          % Excitatory from j to i
                connectivity(newi,newj,2) = nansum(CC(i,j,202:201+20) > repmat(CC_CI(i,j,2),[1,1,length(202:201+20)]));          % Excitatory from i to j
                connectivity(newi,newj,3) = nansum(CC(i,j,201-20:200) < repmat(CC_CI(i,j,1),[1,1,length(201-20:200)]));          % Inhibitory from j to i
                connectivity(newi,newj,4) = nansum(CC(i,j,202:201+20) < repmat(CC_CI(i,j,1),[1,1,length(202:201+20)]));          % Inhibitory from i to j
                
                chanN(newi,newj) = sum(dcp{cellID(i,1,1)}.probe.liveContacts(:));
                fileID(newi,newj) = cellID(i,1,1);
            end
        end
        
        dist = squareform(euclidLoc);
        dist = triu(dist,1);
        nSig3 = triu(nSig2,1);
        nSig3(nSig3 == 0) = NaN;
        for i = 1:4
            connectivity3(:,:,i) = triu(connectivity(:,:,i),1);
        end
        connectivity3(connectivity==0) = NaN;
        chanN = triu(chanN,1);
        fileID = triu(fileID,1);
        
        edges = -0.075:0.15:3.5;
        startID = 1;%find(cellID(:,1,1) > 30,1);
        dist3 = dist(startID:end,startID:end);
        x = nSig3(startID:end,startID:end);
        x = x(:);
        temp = triu(ones(size(dist3)),1);
        temp(temp==0)= NaN;
        dist3 = dist3.*temp;
        dist3 = dist3(:);
        [~,b] = histc(dist3,edges);
        thres = 5:5:50;
        for ti = 1:length(thres)
            for i = 1:length(edges)-1
                probCC(ti,i) = sum(x(b==i) >= thres(ti))/numel(x(b==i));
            end
        end
        
        
        for filei = 1:length(dcp)
            totSig(filei) = nansum(nSig3(fileID==filei));
            totPlus(filei) = sum(nSig3(fileID==filei)>=thres(2));
            chs(filei) = sum(dcp{filei}.probe.liveContacts(:));
        end
        
        figure('Name','Connectivity statistics','Position',[1151 629 1370 700])
        subplot(2,2,[1,3])
        imagesc(nSig3>=thres(2))
        axis square
        xlabel('Neuron #')
        ylabel('Neuron #')
        
        subplot(2,2,2)
        scatter(1:78,totPlus,10*(chs))
        xlabel('Recording day')
        ylabel('Number of pairs')
        
        subplot(2,2,4)
        plot(edges(1:end-1)+(edges(2)-edges(1))/2,probCC')
        xlabel('Euclidean distance')
        ylabel('Average connected pairs')
    end
    
    %% End results file check
end

%% Save results file
if saveResults
    
    saveLocation = ['/home/seth/Projects/DynamicCoherencePhysiology/' dcp{1}.sname ...
        '/initCohDynCohComp'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/neuronTyping' datestr(now,'yyyymmdd')],'-v7.3')
    
end

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

figure
subplot(2,1,1)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(c(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(b(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
end
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,1,:)-a(:,2,2,:)),2)*1000,'Color',initColors(1,:))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,3,:)-a(:,2,2,:)),2)*1000,'Color',initColors(3,:))
plotHorizontal(0);
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Excess spike/s')

subplot(2,1,2)
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,2,:)),2)*1000,'Color',initColors(2,:))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,5,:)),2)*1000,'Color',colors(5,:))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

figure('Name','Grand average','Position',[486 733 1650 420])
minMax = [Inf,0];
for speedi = 1:length(speeds)
    subplot(length(speeds),2,1+(speedi-1)*2)
    for cohi = 1:length(cohs)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speedi,cohi,:)),2)*1000,...
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
    subplot(length(speeds),2,1+(speedi-1)*2)
    ylim(minMax)
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
end

speedi = find(speeds == 10);
subplot(length(speeds),2,2 + (speedi-1)*2)
for seqi = 1:length(sequences)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
    hold on
end
ylim(minMax)
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

gvmrh = figure('Name',['Grand average']);
dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,:),[1,3]))*1000;
for seqi = 1:length(sequences)
    plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,:),[1,3]))*1000;
for seqi = 1:length(sequences)
    plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
end
dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,:),[1,3]))*1000;
for seqi = 1:length(sequences)
    plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
end
if testInitGain
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:),[1,4]))*1000;
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,:),[1,4]))*1000;
            plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
    end
    
end
xlabel('Behavioral gain (unitless)')
ylabel('Spikes/s')

%% Weighted average of MT data
[~,interpolatedWA,mtNeuron_t,~,interpolatedRmt,mtSpref] = mtNeuralAnalysis(...
    'objectFile','allMT_20230419.mat','pltFlg',false,...
    'speeds',[4 8 16 32],'cohs',[10 30 100]);
figure('Name','MT weighted average (interpolated)','Position',[486 733 1650 420])
minMax = [Inf,0];
for speedi = 1:length(speeds)
    subplot(1,length(speeds),speedi)
    for cohi = 1:length(cohs)
        plot(mtNeuron_t,interpolatedWA(:,speedi,cohi),...
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
    subplot(1,length(speeds),speedi)
    ylim(minMax)
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
end

figure('Name',['MT weighted average (interpolated)']);
mtRatesTemp = squeeze(interpolatedWA(mtNeuron_t==450,speeds==10,:));
for seqi = 1:length(sequences)
    plot(gain(seqi,2),mtRatesTemp(cohs == 60),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
mtRatesTemp = squeeze(interpolatedWA(mtNeuron_t==750,speeds==10,:));
mtRatesTemp = [mtRatesTemp(cohs==100) mtRatesTemp(cohs==20) mtRatesTemp(cohs == 100) mtRatesTemp(cohs == 20) mtRatesTemp(cohs == 60)];
for seqi = 1:length(sequences)
    plot(gain(seqi,3),mtRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
end
% mtRatesTemp = squeeze(interpolatedWA(mtNeuron_t==1050,speeds==10,:));
% mtRatesTemp = [mtRatesTemp(cohs==20) mtRatesTemp(cohs==100) mtRatesTemp(cohs == 100) mtRatesTemp(cohs == 20) mtRatesTemp(cohs == 60)];
% for seqi = 1:length(sequences)
%     plot(gain(seqi,4),mtRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
% end
if testInitGain
    for si = 1:length(speeds)
        mtRatesTemp = squeeze(interpolatedWA(mtNeuron_t == 750,si,:));
        plot(squeeze(initGain(si,:,3)),mtRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        hold on
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            mtRatesTemp = squeeze(interpolatedWA(mtNeuron_t == 150,si,:));
            plot(squeeze(initGain(si,:,2)),mtRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
    end
    
end
xlabel('Behavioral gain (unitless)')
ylabel('Spikes/s')

%% Optimize weighed average of integrated MT responses
y = nanmean(Rinit(initCoh.neuron_t>=0 & initCoh.neuron_t<=150,:,:,:),4);
for neuroni = 1:size(interpolatedRmt,4)
    for speedi = 1:length(speeds)
        for cohi = 1:length(cohs)
            tempt = mtNeuron_t(mtNeuron_t >= 0 & mtNeuron_t <= 150);
            xreal = cumsum(interpolatedRmt(mtNeuron_t >= 0 & mtNeuron_t <= 150,speedi,cohi,neuroni));
            xinterp(:,speedi,cohi,neuroni) = interp1(tempt,xreal,0:150);
        end
    end
    X(:,neuroni) = reshape(xinterp(:,:,:,neuroni),[numel(xinterp(:,:,:,neuroni)),1]);
end
liveMT = find(~any(isnan(X),1));
mtWeights = regress(reshape(y,[numel(y),1]),[X(:,liveMT) ones(size(X,1),1)]);
lb = zeros(size(interpolatedRmt,4)+1,1);
ub = [];
mtWeightsnonneg = lsqlin([X(:,liveMT) ones(size(X,1),1)],reshape(y,[numel(y),1]),[],[],[],[],lb,ub);
% interpolatedWAopt = 0 + sum(...
%     repmat(permute(mtWeights(1:end-1),[4,3,2,1]),[sum(mtNeuron_t>=0 & mtNeuron_t<=150),size(interpolatedRmt,2),size(interpolatedRmt,3),1]).*...
%     interpolatedRmt(mtNeuron_t>=0 & mtNeuron_t <=150,:,:,:),4);
iWAopt_t = mtNeuron_t(mtNeuron_t>=0);
interpolatedWAopt = mtWeights(end) + ...
    sum(repmat(permute(mtWeights(1:end-1),[4,3,2,1]),[sum(mtNeuron_t>=0),size(interpolatedRmt,2),size(interpolatedRmt,3),1]) .* ...
    cumsum(interpolatedRmt(mtNeuron_t>=0,:,:,liveMT),1),4);

figure('Name','Optimal weighted sum of temporally integrated MT responses and mean FEF intiation response')
for speedi = 1:length(speeds)
    subplot(3,1,speedi)
    for cohi = 1:length(cohs)
        plot(squeeze(y(:,speedi,:)),'Color',[0 0 0 .2])
        hold on
        plot(mtWeights(end)+squeeze(sum(repmat(permute(mtWeights(1:end-1),[4,3,2,1]),[size(y,1),1,1,size(x,4)]).*xinterp(:,speedi,cohi,liveMT),4)),...
            'Color',[speedColors(speedi,:), cohs(cohi)/100])
    end
    xlabel('Time from motion onset (ms)')
    ylabel('FEF response (spikes/ms)')
end

figure('Name','MT weighted average (interpolated)','Position',[486 733 1650 420])
minMax = [Inf,0];
for speedi = 1:length(speeds)
    subplot(1,length(speeds),speedi)
    for cohi = 1:length(cohs)
        plot(iWAopt_t,interpolatedWAopt(:,speedi,cohi),...
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
    subplot(1,length(speeds),speedi)
    %xlim([0 150])
    ylim(minMax)
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
end

figure('Name',['MT weighted average (interpolated)']);
mtRatesTemp = squeeze(interpolatedWAopt(iWAopt_t==450,speeds==10,:));
for seqi = 1:length(sequences)
    plot(gain(seqi,2),mtRatesTemp(cohs == 60),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
mtRatesTemp = squeeze(interpolatedWAopt(iWAopt_t==750,speeds==10,:));
mtRatesTemp = [mtRatesTemp(cohs==100) mtRatesTemp(cohs==20) mtRatesTemp(cohs == 100) mtRatesTemp(cohs == 20) mtRatesTemp(cohs == 60)];
for seqi = 1:length(sequences)
    plot(gain(seqi,3),mtRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
end
% mtRatesTemp = squeeze(interpolatedWA(mtNeuron_t==1050,speeds==10,:));
% mtRatesTemp = [mtRatesTemp(cohs==20) mtRatesTemp(cohs==100) mtRatesTemp(cohs == 100) mtRatesTemp(cohs == 20) mtRatesTemp(cohs == 60)];
% for seqi = 1:length(sequences)
%     plot(gain(seqi,4),mtRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
% end
if testInitGain
    for si = 1:length(speeds)
        mtRatesTemp = squeeze(interpolatedWAopt(iWAopt_t == 750,si,:));
        plot(squeeze(initGain(si,:,3)),mtRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        hold on
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            mtRatesTemp = squeeze(interpolatedWAopt(iWAopt_t == 150,si,:));
            plot(squeeze(initGain(si,:,2)),mtRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
    end
    
end
xlabel('Behavioral gain (unitless)')
ylabel('Spikes/s')

%% Simple avaraging, by cluster
speedRegression.y = [];
gainRegression.y = [];
if testInitGain
    speedRegression.x = nan(length(speeds)*length(cohs)+3*length(sequences),NumClusters);
    gainRegression.x = nan(length(speeds)*length(cohs)+3*length(sequences),NumClusters);
else
    speedRegression.x = nan(3*length(sequences),NumClusters);
    gainRegression.x = nan(3*length(sequences),NumClusters);
end
for i = 1:NumClusters
    a = Rinit(dynCoh.neuron_t<=1350,:,:,idx == i);
    b = nan(size(Rdyn(:,:,idx == i)));
    c = nan(size(Rdyn(:,:,idx == i)));
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

    figure('Name',['Cluster ' num2str(i)])
    subplot(2,1,1)
    for seqi = 1:5
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(c(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
        hold on
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(b(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
    end
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)-a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(1,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)-a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
    plotHorizontal(0);
    plotVertical([450 750 1050]);
    xlabel('Time from motion onset (ms)')
    ylabel('Excess spike/s')
    
    subplot(2,1,2)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,5,idx == i)),2)*1000,'Color',colors(5,:))
    plotVertical([450 750 1050]);
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    
    figure('Name',['Cluster ' num2str(i)])
    subplot(3,1,1)
    for seqi = 1:5
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(dynCoh.neuron_t<=1350,seqi,idx == i)),2)*1000,...
            'Color',colors(seqi,:))
        hold on
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(dynCoh.neuron_t<=1350,seqi,idx == i)),2)*1000,...
            'Color',colors(seqi,:))
    end
    plotVertical([450 750 1050]);
    xlabel('Time from motion onset (ms)')
    ylabel('Spike/s')
    
    subplot(3,1,2)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)),2)*1000,'Color',initColors(1,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)),2)*1000,'Color',initColors(3,:))
    xlabel('Time from motion onset (ms)')
    ylabel('Spike/s')

    subplot(3,1,3)
    for si = 1:length(speeds)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,si,2,:)),2)*1000,'Color',speedColors(si,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    
    gvrh = figure('Name',['Cluster ' num2str(i)]);
    subplot(1,2,1)
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 450 & eye_tDyn <= 450,:,:),1));
    dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx==i),1)),2)*1000;
    speedRegression.x(1:length(sequences),i) = dynRatesTemp;
    for seqi = 1:length(sequences)
        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
        hold on
    end
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 750 & eye_tDyn <= 750,:,:),1));
    dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx==i),1)),2)*1000;
    speedRegression.x(length(sequences)+1:2*length(sequences),i) = dynRatesTemp;
    for seqi = 1:length(sequences)
        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn(eye_tDyn >= 1050 & eye_tDyn <= 1050,:,:),1));
    dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx==i),1)),2)*1000;
    speedRegression.x(2*length(sequences)+1:3*length(sequences),i) = dynRatesTemp;
    for seqi = 1:length(sequences)
        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    if testInitGain
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
            initRatesTemp = nanmean(squeeze(nanmean(a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:),1)),2)*1000;
            speedRegression.x(3*length(sequences)+(si-1)*length(cohs)+1:3*length(sequences)+si*length(cohs),i) = initRatesTemp;
            plot(eyeSpeedTemp(si,:),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
        if includeEarlyInitCohPertTime
            for si = 1:length(speeds)
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
                initRatesTemp = nanmean(squeeze(nanmean(a(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,:),1)),2)*1000;
                speedRegression.x(3*length(sequences)+3*length(cohs) + (si-1)*length(cohs)+1:3*length(sequences)+3*length(cohs)+si*length(cohs),i) = initRatesTemp;
                plot(eyeSpeedTemp(si,:),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            end
        end
        
    end
    axis square
    set(gca,'TickDir','out')
    xlabel('Eye speed (deg/s)')
    ylabel('Spikes/s')
    
    subplot(1,2,2)
    dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx==i),1)),2)*1000;
    gainRegression.x(1:length(sequences),i) = dynRatesTemp;
    for seqi = 1:length(sequences)
        plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
        hold on
    end
    dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx==i),1)),2)*1000;
    gainRegression.x(length(sequences)+1:2*length(sequences),i) = dynRatesTemp;
    for seqi = 1:length(sequences)
        plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx==i),1)),2)*1000;
    gainRegression.x(2*length(sequences)+1:3*length(sequences),i) = dynRatesTemp;
    for seqi = 1:length(sequences)
        plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    if testInitGain
        for si = 1:length(speeds)
            initRatesTemp = nanmean(squeeze(nanmean(a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:),1)),2)*1000;
            gainRegression.x(3*length(sequences)+(si-1)*length(cohs)+1:3*length(sequences)+si*length(cohs),i) = initRatesTemp;
            plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
        if includeEarlyInitCohPertTime
            for si = 1:length(speeds)
                initRatesTemp = nanmean(squeeze(nanmean(a(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,:),1)),2)*1000;
                gainRegression.x(3*length(sequences)+3*length(cohs) + (si-1)*length(cohs)+1:3*length(sequences)+3*length(cohs)+si*length(cohs),i) = initRatesTemp;
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            end
        end
        
    end
    axis square
    set(gca,'TickDir','out')
    xlabel('Behavioral gain (unitless)')
    ylabel('Spikes/s')
end

%% Gain decoding performance
regressorClusters = 1:NumClusters; %[2,7,8];
gainRegression.x_mean = mean(gainRegression.x,1);
gainRegression.x_std = std(gainRegression.x,[],1);
gainRegression.z = gainRegression.x;
% gainRegression.z = (gainRegression.x - repmat(gainRegression.x_mean,[size(gainRegression.x,1),1])) ./ ...
%     repmat(gainRegression.x_std,[size(gainRegression.x,1),1]);
gainRegression.y(1:length(sequences)) = gain(:,2);
gainRegression.y(length(sequences)+1:2*length(sequences)) = gain(:,3);
gainRegression.y(2*length(sequences)+1:3*length(sequences)) = gain(:,4);
if testInitGain
    for si = 1:length(speeds)
        gainRegression.y(3*length(sequences)+(si-1)*length(cohs)+1:3*length(sequences)+si*length(cohs)) = initGain(si,:,3);
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            gainRegression.y(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+1:3*length(sequences)+3*length(cohs)+si*length(cohs)) = initGain(si,:,2);
        end
    end
end
gainRegression.B = regress(gainRegression.y(:),[gainRegression.z(:,regressorClusters) ones(size(gainRegression.x,1),1)]);
gainRegression.yhat = [gainRegression.z(:,regressorClusters) ones(size(gainRegression.x,1),1)]*gainRegression.B;
gainRegression.Bnonneg = lsqnonneg(gainRegression.z,gainRegression.y');

figure('Name','Gain decoding','Position',[1633 927 888 395])
subplot(1,2,1)
for seqi = 1:length(sequences)
    plot(gainRegression.y(seqi),gainRegression.yhat(seqi),...
        'd','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
for seqi = 1:length(sequences)
    plot(gainRegression.y(seqi+length(sequences)),gainRegression.yhat(seqi+length(sequences)),...
        'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
for seqi = 1:length(sequences)
    plot(gainRegression.y(seqi+2*length(sequences)),gainRegression.yhat(seqi+2*length(sequences)),...
        's','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
if testInitGain
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            plot(gainRegression.y(3*length(sequences)+(si-1)*length(cohs)+ci),...
                gainRegression.yhat(3*length(sequences)+(si-1)*length(cohs)+ci),...
                'o','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
            hold on
        end
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            for ci = 1:length(cohs)
                plot(gainRegression.y(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+ci),...
                    gainRegression.yhat(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+ci),...
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

%% Weighted average of PSTHs according to gain decoding
figure('Name','Weighted average response','Position',[1961 352 577 970])

for i = 1:NumClusters
    a = Rinit(dynCoh.neuron_t<=1350,:,:,idx == i);
    Rinit_weighted_mean(:,:,:,i) = gainRegression.B(i)*(nanmean(a,4)*1000);
    
    Rdyn_weighted_mean(:,:,i) = gainRegression.B(i)*(nanmean(Rdyn(dynCoh.neuron_t<=1350,:,idx == i),3)*1000);
end

for si = 1:length(speeds)
    subplot(length(speeds)+1,1,si)
    for ci = 1:length(cohs)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),sum(Rinit_weighted_mean(:,si,ci,:),4)+gainRegression.B(end),'Color',initColors(ci,:))
        hold on
    end
    axis tight
    ax_wm(si,:) = axis;
    title(['Target speed = ' num2str(speeds(si))] )
end

subplot(length(speeds)+1,1,length(speeds)+1)
for seqi = 1:length(sequences)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),sum(Rdyn_weighted_mean(:,seqi,:),3)+gainRegression.B(end),'Color',colors(seqi,:))
    hold on
end
axis tight
ax_wm(length(si)+1,:) = axis;
title('Dynamic coherence')

for i = 1:length(speeds)+1
    subplot(length(speeds)+1,1,i)
    axis([min(ax_wm(:,1)) max(ax_wm(:,2)) min(ax_wm(:,3)) max(ax_wm(:,4))])
    xlabel('Time from motion onset (ms)')
    ylabel('Estimated gain')
end

figure('Name','Weighted average speed functions','Position',[661 549 1552 281])
subplot(1,4,1)
for ci = 1:length(cohs)
    plot(speeds,sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,:,ci,:),1),4)+gainRegression.B(end),...
        'o-','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
    hold on
end
axis tight
ax_wm2(1,:) = axis;
xlabel('Target speed (deg/s)')
ylabel('Estimated gain')
title('150 ms')

subplot(1,4,2)
for ci = 1:length(cohs)
    plot(speeds,sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,ci,:),1),4)+gainRegression.B(end),...
        'o-','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
    hold on
end
axis tight
ax_wm2(2,:) = axis;
xlabel('Target speed (deg/s)')
ylabel('Estimated gain')
title('450 ms')

subplot(1,4,3)
for ci = 1:length(cohs)
    plot(speeds,sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,ci,:),1),4)+gainRegression.B(end),...
        'o-','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
    hold on
end
axis tight
ax_wm2(3,:) = axis;
xlabel('Target speed (deg/s)')
ylabel('Estimated gain')
title('750 ms')

subplot(1,4,4)
for ci = 1:length(cohs)
    plot(speeds,sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,ci,:),1),4)+gainRegression.B(end),...
        'o-','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:))
    hold on
end
axis tight
ax_wm2(4,:) = axis;
xlabel('Target speed (deg/s)')
ylabel('Estimated gain')
title('1050 ms')

for spi = 1:4
    subplot(1,4,spi)
    axis([min(ax_wm2(:,1)) max(ax_wm2(:,2)) min(ax_wm2(:,3)) max(ax_wm2(:,4))])
end

figure('Name','Weighted average coherence functions','Position',[661 549 1552 281])
subplot(1,4,1)
for si = 1:length(speeds)
    plot(cohs,squeeze(sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,:),1),4)+gainRegression.B(end)),...
        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    hold on
end
axis tight
ax_wm3(1,:) = axis;
xlabel('Coherence')
ylabel('Estimated gain')
title('150 ms')

subplot(1,4,2)
for si = 1:length(speeds)
    plot(cohs,squeeze(sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,si,:,:),1),4)+gainRegression.B(end)),...
        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    hold on
end
axis tight
ax_wm3(2,:) = axis;
xlabel('Coherence')
ylabel('Estimated gain')
title('450 ms')

subplot(1,4,3)
for si = 1:length(speeds)
    plot(cohs,squeeze(sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:),1),4)+gainRegression.B(end)),...
        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    hold on
end
axis tight
ax_wm3(3,:) = axis;
xlabel('Coherence')
ylabel('Estimated gain')
title('750 ms')

subplot(1,4,4)
for si = 1:length(speeds)
    plot(cohs,squeeze(sum(mean(Rinit_weighted_mean(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,si,:,:),1),4)+gainRegression.B(end)),...
        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    hold on
end
axis tight
ax_wm3(4,:) = axis;
xlabel('Coherence')
ylabel('Estimated gain')
title('1050 ms')

for spi = 1:4
    subplot(1,4,spi)
    axis([min(ax_wm3(:,1)) max(ax_wm3(:,2)) min(ax_wm3(:,3)) max(ax_wm3(:,4))])
end

%% Coherence sensitivity vs functional mapping
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

regressionType = 'PCA';
regressorIndex = 1;
switch regressionType
    case 'raw'
        Binittemp = Binit;
        Bdyntemp = Bdyn;
    case 'PCA'
        Binittemp = BinitPCA;
        Bdyntemp = BdynPCA;
end
cax = [NaN NaN];
fh1 = figure('Name','Coherence sensitivity vs functional mapping 750','Position',[627 909 1894 420]);
subplot(1,3,1)
scatter(Y(:,1),Y(:,2),50,mean((Binittemp(regressorIndex+1,:,initCoh.neuron_t>=175 & initCoh.neuron_t<=225)),3),'filled');
cax = caxis;
subplot(1,3,2)
scatter(Y(:,1),Y(:,2),50,mean(Binittemp(regressorIndex+1,:,initCoh.neuron_t>=700 & initCoh.neuron_t<=750),3),'filled');
caxTemp = caxis;
cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];
subplot(1,3,3)
scatter(Y(:,1),Y(:,2),50,mean(Bdyntemp(regressorIndex,:,dynCoh.neuron_t>=700 & dynCoh.neuron_t<=750),3),'filled');
caxTemp = caxis;
cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];



% figure('Name','Coherence \Delta sensitivity vs functional mapping 750','Position',[627 909 1894 420]);
% subplot(1,3,1)
% scatter(Y(:,1),Y(:,2),50,mean((Binit(2,:,initCoh.neuron_t>=175 & initCoh.neuron_t<=225)),3)-mean(Binit(2,:,initCoh.neuron_t>=700 & initCoh.neuron_t<=750),3),'filled');
% caxTemp = caxis;
% cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];
% subplot(1,3,2)
% scatter(Y(:,1),Y(:,2),50,mean((Binit(2,:,initCoh.neuron_t>=175 & initCoh.neuron_t<=225)),3)-mean(Bdyn(1,:,dynCoh.neuron_t>=700 & dynCoh.neuron_t<=750),3),'filled');
% caxTemp = caxis;
% cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];
% subplot(1,3,3)
% scatter(Y(:,1),Y(:,2),50,mean(Binit(2,:,initCoh.neuron_t>=700 & initCoh.neuron_t<=750),3)-mean(Bdyn(1,:,dynCoh.neuron_t>=700 & dynCoh.neuron_t<=750),3),'filled');
% caxTemp = caxis;
% cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];
%
% for i = 1:3
%     subplot(1,3,i)
%     caxis(cax)
%     colorbar
%     xlabel('tSNE_1')
%     ylabel('tSNE_2')
% end



figure('Name','Coherence sensitivity vs functional mapping 1000','Position',[627 909 1894 420])
subplot(1,3,1)
scatter(Y(:,1),Y(:,2),50,mean((Binittemp(regressorIndex+1,:,initCoh.neuron_t>=175 & initCoh.neuron_t<=225)),3),'filled');
caxTemp = caxis;
cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];
subplot(1,3,2)
scatter(Y(:,1),Y(:,2),50,mean(Binittemp(regressorIndex+1,:,initCoh.neuron_t>=950 & initCoh.neuron_t<=1000),3),'filled');
caxTemp = caxis;
cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];
subplot(1,3,3)
scatter(Y(:,1),Y(:,2),50,mean(Bdyntemp(regressorIndex,:,dynCoh.neuron_t>=950 & dynCoh.neuron_t<=1000),3),'filled');
caxTemp = caxis;
cax = [min([cax(1) caxTemp(1)]) max([cax(2) caxTemp(2)])];

for i = 1:3
    subplot(1,3,i)
    caxis(cax)
    colorbar
    xlabel('tSNE_1')
    ylabel('tSNE_2')
end

figure(fh1)
for i = 1:3
    subplot(1,3,i)
    caxis(cax)
    colorbar
    xlabel('tSNE_1')
    ylabel('tSNE_2')
end

% figure('Name','Coherence \Delta sensitivity vs functional mapping 1000','Position',[627 909 1894 420])
% subplot(1,3,1)
% scatter(Y(:,1),Y(:,2),50,mean((Binit(2,:,initCoh.neuron_t>=175 & initCoh.neuron_t<=225)),3)-mean(Binit(2,:,initCoh.neuron_t>=950 & initCoh.neuron_t<=1000),3),'filled');
% cax2(1,:) = caxis;
% subplot(1,3,2)
% scatter(Y(:,1),Y(:,2),50,mean((Binit(2,:,initCoh.neuron_t>=175 & initCoh.neuron_t<=225)),3)-mean(Bdyn(1,:,dynCoh.neuron_t>=950 & dynCoh.neuron_t<=1000),3),'filled');
% cax2(2,:) = caxis;
% subplot(1,3,3)
% scatter(Y(:,1),Y(:,2),50,mean(Binit(2,:,initCoh.neuron_t>=950 & initCoh.neuron_t<=1000),3)-mean(Bdyn(1,:,dynCoh.neuron_t>=950 & dynCoh.neuron_t<=1000),3),'filled');
% cax2(3,:) = caxis;
%
% for i = 1:3
%     subplot(1,3,i)
%     caxis([min(cax2(1,:)) max(cax2(2,:))])
%     colorbar
%     xlabel('tSNE_1')
%     ylabel('tSNE_2')
% end

t0 = 100;
t1 = 750;
t2 = 1050;
alphaLevel = 0.5;
figure('Name','Coherence sensitivity comp','Position',[1237 73 1289 1256])
subplot(3,2,1)
for i = 1:NumClusters
    scatter(Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t1),...
        Bdyntemp(regressorIndex,idx==i,dynCoh.neuron_t==t1),...
        50,'MarkerFaceColor',colorWheel(i,:),'MarkerFaceAlpha',alphaLevel,'MarkerEdgeColor','none')
    hold on
end
axis tight
axTemp(1,:) = axis;
cor = corrcoef(Binittemp(regressorIndex+1,:,initCoh.neuron_t==t1),...
    Bdyntemp(regressorIndex,:,dynCoh.neuron_t==t1));
text(axTemp(1,1),axTemp(1,3),['r^2 = ' num2str(cor(1,2).^2)]);
xlabel(['Init ' num2str(t1)])
ylabel(['Dyn ' num2str(t1)])

subplot(3,2,2)
for i = 1:NumClusters
    scatter(Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t2),...
        Bdyntemp(regressorIndex,idx==i,dynCoh.neuron_t==t2),...
        50,'MarkerFaceColor',colorWheel(i,:),'MarkerFaceAlpha',alphaLevel,'MarkerEdgeColor','none')
    hold on
end
axis tight
axTemp(2,:) = axis;
cor = corrcoef(Binittemp(regressorIndex+1,:,initCoh.neuron_t==t2),...
    Bdyntemp(regressorIndex,:,dynCoh.neuron_t==t2));
text(axTemp(2,1),axTemp(2,3),['r^2 = ' num2str(cor(1,2).^2)]);
xlabel(['Init ' num2str(t2)])
ylabel(['Dyn ' num2str(t2)])

subplot(3,2,3)
for i = 1:NumClusters
    scatter(Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t1),...
        Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t0),...
        50,'MarkerFaceColor',colorWheel(i,:),'MarkerFaceAlpha',alphaLevel,'MarkerEdgeColor','none')
    hold on
end
axis tight
axTemp(3,:) = axis;
cor = corrcoef(Binittemp(regressorIndex+1,:,initCoh.neuron_t==t1),...
    Binittemp(regressorIndex+1,:,initCoh.neuron_t==t0));
text(axTemp(3,1),axTemp(3,3),['r^2 = ' num2str(cor(1,2).^2)]);
xlabel(['Init ' num2str(t1)])
ylabel(['Init ' num2str(t0)])

subplot(3,2,4)
for i = 1:NumClusters
    scatter(Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t2),...
        Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t0),...
        50,'MarkerFaceColor',colorWheel(i,:),'MarkerFaceAlpha',alphaLevel,'MarkerEdgeColor','none')
    hold on
end
axis tight
axTemp(4,:) = axis;
cor = corrcoef(Binittemp(regressorIndex+1,:,initCoh.neuron_t==t2),...
    Binittemp(regressorIndex+1,:,initCoh.neuron_t==t0));
text(axTemp(4,1),axTemp(4,3),['r^2 = ' num2str(cor(1,2).^2)]);
xlabel(['Init ' num2str(t2)])
ylabel(['Init ' num2str(t0)])

subplot(3,2,5)
for i = 1:NumClusters
    scatter(Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t1),...
        Binittemp(regressorIndex+1,idx==i,initCoh.neuron_t==t2),...
        50,'MarkerFaceColor',colorWheel(i,:),'MarkerFaceAlpha',alphaLevel,'MarkerEdgeColor','none')
    hold on
end
axis tight
axTemp(5,:) = axis;
cor = corrcoef(Binittemp(regressorIndex+1,:,initCoh.neuron_t==t1),...
    Binittemp(regressorIndex+1,:,initCoh.neuron_t==t2));
text(axTemp(5,1),axTemp(5,3),['r^2 = ' num2str(cor(1,2).^2)]);
xlabel(['Init ' num2str(t1)])
ylabel(['Init ' num2str(t2)])

subplot(3,2,6)
for i = 1:NumClusters
    scatter(Bdyntemp(regressorIndex,idx==i,dynCoh.neuron_t==t1),...
        Bdyntemp(regressorIndex,idx==i,dynCoh.neuron_t==t2),...
        50,'MarkerFaceColor',colorWheel(i,:),'MarkerFaceAlpha',alphaLevel,'MarkerEdgeColor','none')
    hold on
end
axis tight
axTemp(6,:) = axis;
cor = corrcoef(Bdyntemp(regressorIndex,:,dynCoh.neuron_t==t1),...
    Bdyntemp(regressorIndex,:,dynCoh.neuron_t==t2));
text(axTemp(6,1),axTemp(6,3),['r^2 = ' num2str(cor(1,2).^2)]);
xlabel(['Dyn ' num2str(t1)])
ylabel(['Dyn ' num2str(t2)])

for i = 1:6
    subplot(3,2,i)
    axis([min(axTemp(:)) max(axTemp(:)) min(axTemp(:)) max(axTemp(:))])
    plotUnity;
    axis square
    plotVertical(0);
    plotHorizontal(0);
end

%% Make movie
movieType = 'PCA';
movieDir = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/';
switch movieType
    case 'raw'
        for ni = 1:size(Binit,2)
            cdat(ni,:) = smooth(Binit(2,ni,:),5);
        end
    case 'PCA'
        for ni = 1:size(BinitPCA,2)
            cdat(ni,:) = smooth(BinitPCA(2,ni,:),5);
        end

end
cax = [min(cdat,[],[1,2]) max(cdat,[],[1,2])];
hm = figure('Name','Sensitivty movie','Position',[1961 339 560 990]);
ti = 1;
subplot(2,1,1)
imagesc(initCoh.neuron_t,1:212,cdat)
hold on
hv = plotVertical(initCoh.neuron_t(ti));
xlabel('Time from motion onset (ms)')
ylabel('Neuron #')
set(hv,'XData',[initCoh.neuron_t(ti) initCoh.neuron_t(ti)])
subplot(2,1,2)
hc = scatter(Y(:,1),Y(:,2),50,cdat(:,ti),'filled');
caxis(cax)
axis square
xlabel('tSNE_1')
ylabel('tSNE_2')
for ti = 1:size(BinitPCA,3)
    set(hv,'XData',[initCoh.neuron_t(ti) initCoh.neuron_t(ti)])
    set(hc,'cdata',cdat(:,ti))
    caxis(cax)
    frame = getframe(hm);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ti == 1
        imwrite(imind,cm,[movieDir movieType 'Animation.gif'],'gif', 'Loopcount',inf,'DelayTime',1/100);
    else
        imwrite(imind,cm,[movieDir movieType 'Animation.gif'],'gif','WriteMode','append','DelayTime',1/100);
    end
end


%% PCA plotting

figure
for ri = 1:reconstructionN
    for si = 1:size(RlowD,2)
        subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
        for ci = 1:size(RlowD,3)
            plot(initCoh.neuron_t,RlowD(:,si,ci,ri),'Color',initColors(ci,:))
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
% subplot(2,2,1)
% for cohi = 1:3
%     plot(initCoh.neuron_t,pInit(:,2,cohi,2),'-','Color',initColors(cohi,:))
%     hold on
% end
% for seqi = 1:5
%     plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,2),'-','Color',colors(seqi,:))
%     hold on
% end
% axis tight
% ax(1,:) = axis;
% subplot(2,2,2)
% for cohi = 1:3
%     plot(initCoh.neuron_t,pInit(:,2,cohi,2)-pInit(:,2,2,2),'-','Color',initColors(cohi,:))
%     hold on
% end
% for seqi = 1:5
%     plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,2)-pDynCross(:,5,2),'-','Color',colors(seqi,:))
%     hold on
% end
% axis tight
%
% subplot(2,2,3)
% for cohi = 1:3
%     plot(initCoh.neuron_t,pInit(:,2,cohi,1),'-','Color',initColors(cohi,:))
%     hold on
% end
% for seqi = 1:5
%     plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,1),'-','Color',colors(seqi,:))
%     hold on
% end
% axis tight
%
% subplot(2,2,4)
% for cohi = 1:3
%     plot(initCoh.neuron_t,pInit(:,2,cohi,1)-pInit(:,2,2,1),'-','Color',initColors(cohi,:))
%     hold on
% end
% for seqi = 1:5
%     plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,1)-pDynCross(:,5,1),'-','Color',colors(seqi,:))
%     hold on
% end
% axis tight

subplot(3,3,1)
for speedi = 1:3
    for cohi = 1:3
        plot(initCoh.neuron_t,pInit(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
axis tight
ax(1,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Speed related activity (a.u.)')

subplot(3,3,2)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,1),'-','Color',colors(seqi,:))
    hold on
end
axis tight
ax(2,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Speed related activity (a.u.)')

subplot(3,3,6)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDyn(:,seqi,1),'-','Color',colors(seqi,:))
    hold on
end
axis tight
xlabel('Time from motion onset (ms)')
ylabel('Coherence related activity (a.u.)')
% ax(3,:) = axis;

for subploti = 1:2
    subplot(3,3,subploti)
    axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    if subploti == 1
        plotVertical(75);
    elseif subploti == 3
        plotVertical(750);
    end
end

subplot(3,3,4)
for cohi = 1:3
    for speedi = 1:length(speeds)
        plot(initCoh.neuron_t,pInit(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
axis tight
ax(1,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Coherence related activity (a.u.)')

subplot(3,3,5)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,2),'-','Color',colors(seqi,:))
    hold on
end
axis tight
ax(2,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Coherence related activity (a.u.)')

for subploti = 1:2
    subplot(3,3,subploti+3)
    axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    if subploti == 1
        plotVertical(75);
    elseif subploti == 3
        plotVertical(750);
    end
end

subplot(3,3,7)
for cohi = 1:3
    for speedi = 1:length(speeds)
        plot(initCoh.neuron_t,pInit(:,speedi,cohi,3),'-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
axis tight
ax(1,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Gain related activity (a.u.)')

subplot(3,3,8)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDynCross(:,seqi,3),'-','Color',colors(seqi,:))
    hold on
end
axis tight
ax(2,:) = axis;
xlabel('Time from motion onset (ms)')
ylabel('Gain related activity (a.u.)')

for subploti = 1:2
    subplot(3,3,subploti+6)
    axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    if subploti == 1
        plotVertical(75);
    elseif subploti == 3
        plotVertical(750);
    end
end

subplot(3,3,9)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pDyn(:,seqi,2),'-','Color',colors(seqi,:))
    hold on
end
axis tight
xlabel('Time from motion onset (ms)')
ylabel('Gain related activity (a.u.)')

%% Targeted dimension subspace
figure('Name','Targeted subspace')
subplot(1,3,1)
initTaccept = initCoh.neuron_t>=-100 & initCoh.neuron_t<=1350;
dynTaccept = dynCoh.neuron_t>=-100 & dynCoh.neuron_t<=1350;
for speedi = 1:3
    for cohi = 1:3
        plot(...
            pInit(initTaccept,speedi,cohi,2),pInit(initTaccept,speedi,cohi,3),...
            '-','Color',[speedColors(speedi,:) cohs(cohi)/100])
        hold on
    end
end
grid on
xlabel('Speed related activity')
ylabel('Gain related activity')

subplot(1,3,2)
for seqi = 1:5
    plot(...
        pDynCross(dynTaccept,seqi,2),pDynCross(dynTaccept,seqi,3),...
        '-','Color',colors(seqi,:))
    hold on
end
grid on
xlabel('Speed related activity')
ylabel('Gain related activity')

subplot(1,3,3)
for seqi = 1:5
    plot(pDyn(dynTaccept,seqi,1),pDyn(dynTaccept,seqi,2),'-','Color',colors(seqi,:))
    hold on
end
xlabel('Coherence related activity')
ylabel('Gain related activity')

%% Targeted dimension activity vs gain
dimNames = {'Speed','Coherence','Gain','Offset'};
gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
for targetedDim = 1:3
    tempData = [];
    subplot(3,1,targetedDim)
    dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
        hold on
    end
    tempData = [tempData; gain(:,2),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    tempData = [tempData; gain(:,3),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
       
    if testInitGain
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pInit(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
        end
        if includeEarlyInitCohPertTime
            for si = 1:length(speeds)
                initRatesTemp = squeeze(nanmean(pInit(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
        end
        
    end
    gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
    axis square
    ax = axis;
    text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
    xlabel('Behavioral gain (unitless)')
    ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
end


gvth2 = figure('Name',['Behavioral gain vs activity on targeted dimension (dynCoh)'],'Position',[1956 59 570 1263]);
for targetedDim = 2:3
    tempData = [];
    subplot(3,1,targetedDim)
    dynRatesTemp = squeeze(nanmean(pDyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim-1),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
        hold on
    end
    tempData = [tempData; gain(:,2),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pDyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim-1),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    tempData = [tempData; gain(:,3),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pDyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim-1),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    tempData = [tempData; gain(:,4),dynRatesTemp(:)];
    
    gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
    axis square
    ax = axis;
    text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
    xlabel('Behavioral gain (unitless)')
    ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
end

%% Strength of inputs in targeted dimension
figure('Name','Motion input strength','Position',[1440 757 1086 572])
subplot(2,2,1)
plot(initCoh.neuron_t,betaInit(1,:,1))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),betaDyn(1,:,1))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Speed related sensitivity (a.u.)')

subplot(2,2,2)
plot(initCoh.neuron_t,betaInit(2,:,1))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),betaDyn(2,:,1))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Speed related input (a.u.)')

subplot(2,2,3)
plot(initCoh.neuron_t,betaInit(1,:,2))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),betaDyn(1,:,2))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Coherence related sensitivity (a.u.)')

subplot(2,2,4)
plot(initCoh.neuron_t,betaInit(2,:,2))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),betaDyn(2,:,2))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Coherence related input (a.u.)')

%%
figure('Name','State vs coherence sensitivity','Position',[1178 909 1343 420])
subplot(1,2,1)
indvec = initCoh.neuron_t >= 550 & initCoh.neuron_t <= 1150;
tempR = RlowD(indvec,2,2,:);
tempBeta(:,1) = betaInit(1,indvec,2);
indvec = dynCoh.neuron_t >= 550 & dynCoh.neuron_t <= 1150;
tempBeta(:,2) = betaDyn(1,indvec,2);
tempBeta = (tempBeta-min(tempBeta(:)))/(max(tempBeta(:))-min(tempBeta(:)));
pcdims = [1,2];
for ti = 1:size(tempR,1)
    plot(tempR(ti,1,1,pcdims(1)),tempR(ti,1,1,pcdims(2)),'o','Color',[0 0 tempBeta(ti,1)],...
        'MarkerFaceColor',[0 0 tempBeta(ti,1)])
    hold on
end
tempR = RDynlowD(indvec,5,:);
for ti = 1:size(tempR,1)
    plot(tempR(ti,1,pcdims(1)),tempR(ti,1,pcdims(2)),'o','Color',[tempBeta(ti,2) 0 0],...
        'MarkerFaceColor',[tempBeta(ti,2) 0 0])
    hold on
end
% scatter(RlowD(indvec,2,2,1),RlowD(indvec,2,2,2),15,betaInit(1,indvec,2),'filled')
% hold on
% indvec = dynCoh.neuron_t >= 550 & dynCoh.neuron_t <= 1150;
% scatter(RDynlowD(indvec,5,1),RDynlowD(indvec,5,2),15,betaDyn(1,indvec,2),'filled')
indvec = initCoh.neuron_t == 750;
indvec2 = dynCoh.neuron_t == 750;
plot([RlowD(indvec,2,2,pcdims(1)) RDynlowD(indvec,5,pcdims(1))],[RlowD(indvec,2,2,pcdims(2)) RDynlowD(indvec,5,pcdims(2))],'ko-')
indvec = initCoh.neuron_t == 1050;
indvec2 = dynCoh.neuron_t == 1050;
plot([RlowD(indvec,2,2,pcdims(1)) RDynlowD(indvec,5,pcdims(1))],[RlowD(indvec,2,2,pcdims(2)) RDynlowD(indvec,5,pcdims(2))],'ko-')
xlabel('PC 1')
ylabel('PC 2')
axis square
%colorbar

subplot(1,2,2)
tempBeta2 = nan(size(betaInit,2),2);
tempBeta2(:,1) = betaInit(1,:,2);
tempBeta2(1:size(betaDyn,2),2) = betaDyn(1,:,2);
tempBeta2 = (tempBeta2-min(tempBeta2(:)))/(max(tempBeta2(:))-min(tempBeta2(:)));
for ti = 1:2:size(RlowD,1)
    plot(RlowD(ti,2,2,pcdims(1)),RlowD(ti,2,2,pcdims(2)),'o','Color',[0 0 tempBeta2(ti,1)],...
        'MarkerFaceColor',[0 0 tempBeta2(ti,1)])
    hold on
end
for ti = 1:2:size(RDynlowD,1)
    if ~isnan(tempBeta2(ti,2))
        plot(RDynlowD(ti,5,pcdims(1)),RDynlowD(ti,5,pcdims(2)),'o','Color',[tempBeta2(ti,2) 0 0],...
            'MarkerFaceColor',[tempBeta2(ti,2) 0 0])
    end
end
% scatter(RlowD(:,2,2,1),RlowD(:,2,2,2),15,betaInit(1,:,2),'filled')
% hold on
% scatter(RDynlowD(:,5,1),RDynlowD(:,5,2),15,betaDyn(1,:,2),'filled')
indvec = initCoh.neuron_t == 750;
indvec2 = dynCoh.neuron_t == 750;
plot([RlowD(indvec,2,2,pcdims(1)) RDynlowD(indvec,5,pcdims(1))],[RlowD(indvec,2,2,pcdims(2)) RDynlowD(indvec,5,pcdims(2))],'ko-')
indvec = initCoh.neuron_t == 1050;
indvec2 = dynCoh.neuron_t == 1050;
plot([RlowD(indvec,2,2,pcdims(1)) RDynlowD(indvec,5,pcdims(1))],[RlowD(indvec,2,2,pcdims(2)) RDynlowD(indvec,5,pcdims(2))],'ko-')
xlabel('PC 1')
ylabel('PC 2')
axis square
% colorbar
% indvec = initCoh.neuron_t >= 450 & initCoh.neuron_t <= 1050;
% scatter(RlowD(indvec,2,2,3),RlowD(indvec,2,2,2),15,betaInit(1,indvec,2))
% hold on
% indvec = dynCoh.neuron_t >= 450 & dynCoh.neuron_t <= 1050;
% scatter(RDynlowD(indvec,5,3),RDynlowD(indvec,5,2),15,betaDyn(1,indvec,2))
% indvec = initCoh.neuron_t == 740;
% indvec2 = dynCoh.neuron_t == 740;
% plot([RlowD(indvec,2,2,3) RDynlowD(indvec,5,3)],[RlowD(indvec,2,2,2) RDynlowD(indvec,5,2)],'ko-')
% indvec = initCoh.neuron_t == 1040;
% indvec2 = dynCoh.neuron_t == 1040;
% plot([RlowD(indvec,2,2,3) RDynlowD(indvec,5,3)],[RlowD(indvec,2,2,2) RDynlowD(indvec,5,2)],'ko-')

%% Plot activity along orthonormal subspace defined by targeted optimized across time points
figure('Name','Activity along targeted dimensions found by optimizing over time',...
    'Position',[73 164 1616 1072])
dimensionLabels = {'Speed','Coherence','Gain'};
for di = 1:3
    optH(1+(di-1)*3) = subplot(3,3,1+(di-1)*3);
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            plot(initCoh.neuron_t,pBopt(:,si,ci,di),'-','Color',[speedColors(si,:) cohs(ci)/100])
            hold on
        end
    end
    axis tight
    xlabel('Time from motion onset (ms)')
    ylabel('Activity along targeted dimension')
    title(dimensionLabels{di})
    
    optH(2+(di-1)*3) = subplot(3,3,2+(di-1)*3);
    for seqi = 1:length(sequences)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDynCross(:,seqi,di),'-','Color',colors(seqi,:))
        hold on
    end
    axis tight
    xlabel('Time from motion onset (ms)')
    ylabel('Activity along targeted dimension')
    title(dimensionLabels{di})
end

optH(6) = subplot(3,3,6);
for seqi = 1:length(sequences)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDyn(:,seqi,1),'-','Color',colors(seqi,:))
    hold on
end
    axis tight
xlabel('Time from motion onset (ms)')
ylabel('Activity along targeted dimension')
title(dimensionLabels{2})

optH(9) = subplot(3,3,9);
for seqi = 1:length(sequences)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),pOptDyn(:,seqi,2),'-','Color',colors(seqi,:))
    hold on
end
    axis tight
xlabel('Time from motion onset (ms)')
ylabel('Activity along targeted dimension')
title(dimensionLabels{3})

linkaxes(optH,'xy')

%% Plot behavioral gain vs activity along targeted dimensions
dimNames = {'Speed','Coherence','Gain','Offset'};
gvthOpt = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
for targetedDim = 1:3
    tempData = [];
    subplot(3,1,targetedDim)
    dynRatesTemp = squeeze(nanmean(pOptDynCross(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
        hold on
    end
    tempData = [tempData; gain(:,2),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pOptDynCross(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    tempData = [tempData; gain(:,3),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pOptDynCross(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
       
    if testInitGain
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(pBopt(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,targetedDim),1));
            plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
        end
        if includeEarlyInitCohPertTime
            for si = 1:length(speeds)
                initRatesTemp = squeeze(nanmean(pBopt(dynCoh.neuron_t >= 100 & dynCoh.neuron_t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
        end
        
    end
    gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
    axis square
    ax = axis;
    text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
    xlabel('Behavioral gain (unitless)')
    ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
end


gvthOpt2 = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (dynCoh)'],'Position',[1956 59 570 1263]);
for targetedDim = 2:3
    tempData = [];
    subplot(3,1,targetedDim)
    dynRatesTemp = squeeze(nanmean(pOptDyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim-1),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
        hold on
    end
    tempData = [tempData; gain(:,2),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pOptDyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim-1),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    tempData = [tempData; gain(:,3),dynRatesTemp(:)];
    dynRatesTemp = squeeze(nanmean(pOptDyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim-1),1));
    for seqi = 1:length(sequences)
        plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    end
    tempData = [tempData; gain(:,4),dynRatesTemp(:)];
    
    gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
    axis square
    ax = axis;
    text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
    xlabel('Behavioral gain (unitless)')
    ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
end

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
   figure('Name',['Unit ' dcp{exIndex(1)}.datapath(end-11:end) '.' num2str(exIndex(2))],...
       'Position',[1770 1103 751 226])

   subplot(1,2,1)
   for ci = 1:3
       plot(initCoh.neuron_t,Rinit(:,2,ci,listIndex)*1000,'Color',initColors(ci,:));
       hold on
   end
   xlabel('Time from motion onset (ms)')
   ylabel('Spikes/s')
   axis tight
   ax(1,:) = axis;

   subplot(1,2,2)
   for seqi = 1:5
       plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,seqi,listIndex)*1000,'Color',colors(seqi,:));
       hold on
   end
   xlabel('Time from motion onset (ms)')
   ylabel('Spikes/s')
   axis tight
   ax(2,:) = axis;

   newax = [min(ax(:,1)),max(ax(:,2)),min(ax(:,3)),max(ax(:,4))];
   subplot(1,2,1)
   axis(newax)
   subplot(1,2,2)
   axis(newax)
   text((newax(2)-newax(1))*.5+newax(1),(newax(4)-newax(3))*.1+newax(3),...
       [dcp{exIndex(1)}.datapath(end-11:end) '.' num2str(exIndex(2))])
end

%%
figure('Name','F neuron 2','Position',[1770 1103 751 226])
exIndex = [67, 186];
cellID2 = squeeze(cellID(:,1,1:2));
listIndex = find(ismember(cellID2, exIndex, 'rows'));

subplot(1,2,1)
for ci = 1:3
    plot(initCoh.neuron_t,Rinit(:,2,ci,listIndex)*1000,'Color',initColors(ci,:));
    hold on
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')
axis tight
ax(1,:) = axis;

subplot(1,2,2)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,seqi,listIndex)*1000,'Color',colors(seqi,:));
    hold on
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')
axis tight
ax(2,:) = axis;

subplot(1,2,1)
axis([min(ax(:,1)),max(ax(:,2)),min(ax(:,3)),max(ax(:,4))])
subplot(1,2,2)
axis([min(ax(:,1)),max(ax(:,2)),min(ax(:,3)),max(ax(:,4))])

figure('Name','Early vs late neuron','Position',[1770 1103 751 226])
exIndex = [69, 81];
cellID2 = squeeze(cellID(:,1,1:2));
listIndex = find(ismember(cellID2, exIndex, 'rows'));

subplot(1,2,1)
for ci = 1:3
    plot(initCoh.neuron_t,Rinit(:,2,ci,listIndex)*1000,'Color',initColors(ci,:));
    hold on
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')
axis tight
ax(1,:) = axis;

subplot(1,2,2)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,seqi,listIndex)*1000,'Color',colors(seqi,:));
    hold on
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')
axis tight
ax(2,:) = axis;

subplot(1,2,1)
axis([min(ax(:,1)),max(ax(:,2)),min(ax(:,3)),max(ax(:,4))])
subplot(1,2,2)
axis([min(ax(:,1)),max(ax(:,2)),min(ax(:,3)),max(ax(:,4))])

%% Visualize dimensionality reduction
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
        plot(initCoh.neuron_t,...
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
        plot(initCoh.neuron_t,...
            nanmean(squeeze(Rinit(:,3,3,idx == i))./...
            repmat(max(squeeze(Rinit(:,3,3,idx==i)),[],1),[size(Rinit,1) 1]),2),...
            'Color',colorWheel(i,:))
        hold on
    end

    subplot(4,2,[7,8])
    for i = 1:NumClusters
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
            nanmean(squeeze(Rdyn(:,5,idx == i))./...
            repmat(max(squeeze(Rdyn(:,5,idx==i)),[],1),[size(Rdyn,1) 1]),2),...
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
            plot(initCoh.neuron_t,...
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

figure('Name','Coherence effect, dynCoh','Position',[1956 196 570 1133])
for i = 1:NumClusters
    subplot(NumClusters,1,i)
    for seqi = 1:size(Rdyn,2)
        tempR = nanmean(squeeze(Rdyn(:,seqi,idx == i))./...
            repmat(max(squeeze(Rdyn(:,seqi,idx==i)),[],1),[size(Rdyn,1) 1]),2);
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
            tempR,...
            'Color',colors(seqi,:))
        hold on
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
    plot(initCoh.neuron_t(initCoh.neuron_t<=1350),...
        nanmean((squeeze(Rinit(initCoh.neuron_t<=1350,3,3,idx == i))-nanmean(squeeze(Rinit(initCoh.neuron_t<=1350,2,3,idx == i)),1))./...
        repmat(nanstd(squeeze(Rinit(initCoh.neuron_t<=1350,2,3,idx==i)),[],1),[size(Rinit(initCoh.neuron_t<=1350,:,:,:),1) 1]),2),...
        'Color',colorWheel(i,:))
    xlabel('Time from motion onset (ms)')
    ylabel('z-score')
end

%% Plot zscores by 1D projection of functional topology
thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
[~,thetaSort] = sort(thetas);

figure
subplot(1,2,1)
imagesc(initCoh.neuron_t,1:size(zInit,1),squeeze(zInit(:,2,2,thetaSort))')
xlabel('Time from motion onset (ms)')
ylabel('Sorted neuron #')

subplot(1,2,2)
imagesc(dynCoh.neuron_t,1:size(zDyn,1),squeeze(zDyn(:,5,thetaSort))')
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
%         zR = (squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))-mean(squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))))/std(squeeze(Rinit(:,3,3,neighbors(exUnit,ni))));
%         plot(initCoh.neuron_t,zR,'Color',colors(ni,:))
        zR = (squeeze(Rdyn(:,5,neighbors(exUnit,ni)))-mean(squeeze(Rdyn(:,5,neighbors(exUnit,ni)))))/std(squeeze(Rdyn(:,5,neighbors(exUnit,ni))));
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),zR,'Color',colors(ni,:))
        hold on
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

%% Connection inferred from CC

if calcCC
    % Measure connectivity
    thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
    [xtemp,ytemp] = generalEllipse(1,1,'theta',thetas);

    figure
    circleProject = false;
    plotPairDetail = false;
    thresVal = 1;
    piOr2pi = '2pi';

    if circleProject
        scatter(xtemp,ytemp,50,colorWheel(idx,:),'filled');
        hold on

        for uniti = 1:size(nSig3,1)
            for unitj = 1:size(nSig3,2)
                if ~isnan(nSig3(uniti,unitj)) & nSig3(uniti,unitj) >= thresVal
                    plot([xtemp(uniti) xtemp(unitj)],[ytemp(uniti) ytemp(unitj)],'k-',...
                        'lineWidth',10*(nSig3(uniti,unitj)/max(nSig3(:)))^2)
                end
            end
        end

        xlabel('circle X')
        ylabel('circle Y')
    else
        axInd = 1;
        for i = 1:4
            subplot(1,2,axInd)
            if i == 1 | i == 3
                scatter(Y(:,1),Y(:,2),50,colorWheel(idx,:),'filled');
            else
                axInd = axInd+1;
            end
            hold on

            for uniti = 1:size(nSig3,1)
                for unitj = uniti:size(nSig3,2)
                    if ~isnan(connectivity3(uniti,unitj,i)) & connectivity3(uniti,unitj,i) >= thresVal
                        if i == 1 | i == 3
                            quiver(Y(unitj,1),Y(unitj,2),Y(uniti,1)-Y(unitj,1),Y(uniti,2)-Y(unitj,2),0,'k',...
                                'lineWidth',5*(connectivity(uniti,unitj,i)/max(connectivity3(:))))
                        else
                            quiver(Y(uniti,1),Y(uniti,2),Y(unitj,1)-Y(uniti,1),Y(unitj,2)-Y(uniti,2),0,'k',...
                                'lineWidth',5*(connectivity(uniti,unitj,i)/max(connectivity3(:))))
                        end
                        %                         plot([Y(uniti,1) Y(unitj,1)],[Y(uniti,2) Y(unitj,2)],'k-',...
                        %                             'lineWidth',10*(nSig3(uniti,unitj)/max(nSig3(:)))^2)
                    end
                end
            end

            xlabel('tSNE1')
            ylabel('tSNE2')
            axis square
        end
    end

    figure
    axInd = 1;
    for i = 1:4
        subplot(1,2,axInd)
        if i == 1 | i == 3

        else
            axInd = axInd+1;
        end
        hold on
        tempXYS = [];
        for uniti = 1:size(nSig3,1)
            for unitj = uniti:size(nSig3,2)
                if ~isnan(connectivity3(uniti,unitj,i)) & connectivity3(uniti,unitj,i) >= thresVal
                    if i == 1 | i == 3
                        tempXYS = [tempXYS; thetas(unitj),thetas(uniti),connectivity(uniti,unitj,i)/max(connectivity(:)),idx(uniti)];
                    else
                        tempXYS = [tempXYS; thetas(uniti),thetas(unitj),connectivity(uniti,unitj,i)/max(connectivity(:)),idx(unitj)];
                    end
                end
            end
        end
        switch piOr2pi
            case 'pi'
                scatter((tempXYS(:,1)),(tempXYS(:,2)),100*tempXYS(:,3).^2,colorWheel(tempXYS(:,4),:),'filled')
            case '2pi'
                scatter(wrapTo2Pi(tempXYS(:,1)),wrapTo2Pi(tempXYS(:,2)),100*tempXYS(:,3).^2,colorWheel(tempXYS(:,4),:),'filled')
        end
    end
    for subploti = 1:2
        subplot(1,2,subploti)
        switch piOr2pi
            case 'pi'
                axis([-pi pi -pi pi])
                plotHorizontal([-pi/2 0 pi/2]);
                plotVertical([-pi/2 0 pi/2]);

            case '2pi'
                axis([0 2*pi 0 2*pi])
                plotHorizontal([pi/2 pi 3*pi/2]);
                plotVertical([pi/2 pi 3*pi/2]);
        end
        axis square

        xlabel('\theta_{pre}')
        ylabel('\theta_{post}')
    end

    if plotPairDetail
        for uniti = 1:size(nSig3,1)
            for unitj = 1:size(nSig3,2)
                for i = 1:4
                    if ~isnan(connectivity3(uniti,unitj,i)) & connectivity3(uniti,unitj,i) >= thresVal
                        figure('Name','Pairs in detail','Position',[212 337 1430 750])
                        subplot(3,2,[1,3,5])
                        scatter(Y(:,1),Y(:,2),50,colorWheel(idx,:),'filled');
                        hold on
                        if i == 1 | i == 3
                            quiver(Y(unitj,1),Y(unitj,2),Y(uniti,1)-Y(unitj,1),Y(uniti,2)-Y(unitj,2),0,'k',...
                                'lineWidth',5*(connectivity(uniti,unitj,i)/max(connectivity3(:))))
                        else
                            quiver(Y(uniti,1),Y(uniti,2),Y(unitj,1)-Y(uniti,1),Y(unitj,2)-Y(uniti,2),0,'k',...
                                'lineWidth',5*(connectivity(uniti,unitj,i)/max(connectivity3(:))))
                        end
                        plot(Y(uniti,1),Y(uniti,2),'kx','MarkerSize',10)
                        plot(Y(unitj,1),Y(unitj,2),'ks','MarkerSize',10)
                        axis square
                        xlabel('tSNE1')
                        ylabel('tSNE2')

                        subplot(3,2,2)
                        j = find(cellID(uniti,:,3)==cellID(unitj,1,2));
                        plot(-200:200,squeeze(CC(uniti,j,:)))
                        hold on
                        plotHorizontal(CC_CI(uniti,j,:));
                        plotVertical(0);
                        lims = axis;
                        text(0.8*(lims(2)-lims(1))+lims(1),lims(4),[num2str(nSig3(uniti,unitj))])
                        xlabel(['Time from spike neuron '  num2str(cellID(uniti,1,2)) ,' file ' num2str(cellID(uniti,1,1))])
                        ylabel(['Excess spikes neuron ' num2str(cellID(uniti,j,3)) ', file ' num2str(cellID(uniti,1,1))])

                        subplot(3,2,4)
                        plot(initCoh.neuron_t,Rinit(:,2,2,uniti)*1000)
                        hold on
                        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,5,uniti)*1000)
                        xlabel('Time from motion onset (ms)')
                        ylabel('Spikes/s')

                        subplot(3,2,6)
                        cell2indx = find(cellID(:,1,1) == cellID(uniti,1,1) & cellID(:,1,2) == cellID(uniti,j,3));
                        plot(initCoh.neuron_t,Rinit(:,2,2,cell2indx)*1000)
                        hold on
                        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,5,cell2indx)*1000)
                        xlabel('Time from motion onset (ms)')
                        ylabel('Spikes/s')

                        input('Enter to continue')
                    end
                end
            end
        end
    end


end
