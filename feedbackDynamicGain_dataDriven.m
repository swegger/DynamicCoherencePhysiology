function feedbackDynamicGain_dataDriven(varargin)
%%
%
%
%
%%

%% Defaults


%% Parse inputs
%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[5,10,20])
addParameter(Parser,'cohs',[20, 60, 100])
addParameter(Parser,'dir',[0 180])
addParameter(Parser,'initWin',[150 200])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'resultsFile','none')
addParameter(Parser,'testInitGain',true)
addParameter(Parser,'includeEarlyInitCohPertTime',false)
addParameter(Parser,'sequences',[1; 2; 3; 4; 5])
addParameter(Parser,'perturbations',[0; 4; 6; 8])
addParameter(Parser,'NumClusters',8)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'Gfb',[0.995 0.997 0.999])

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
dir = Parser.Results.dir;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
resultsFile = Parser.Results.resultsFile;
testInitGain = Parser.Results.testInitGain;
includeEarlyInitCohPertTime = Parser.Results.includeEarlyInitCohPertTime;
sequences = Parser.Results.sequences;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
Gfb = Parser.Results.Gfb;


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

%% Results file check start
if ~fileExist
    
    disp('Neural analysis loop...')
    
    %% Get mean and covariance of each unit
    passCutoff = nan(1000,1);
    Rinit = nan(1701,3,3,1000);
    Rdyn = nan(1701,5,1000);
    locations = nan(1000,3);
    cellID = nan(1000,100,3);
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
    dynCohT = [60*ones(sum(dynCoh.neuron_t<0),5);dynCohT];
    % Linear model
    Binit = nan(4,size(zInit,4),size(zInit,1));
    Bdyn = nan(3,size(zInit,4),size(zDyn,1));
    lagT = 100;
    [Cohs, Spds] = meshgrid(cohs,speeds);
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
    
    
    %% End results file check
end

%% Convert FEF gain signal to behavioral gain

targetedDim = 3;
tempData = [];
dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,targetedDim),1));
tempData = [tempData; gain(:,2),dynRatesTemp(:)];
dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,targetedDim),1));
tempData = [tempData; gain(:,3),dynRatesTemp(:)];
dynRatesTemp = squeeze(nanmean(pDynCross(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,targetedDim),1));
tempData = [tempData; gain(:,4),dynRatesTemp(:)];

for si = 1:length(speeds)
    initRatesTemp = squeeze(nanmean(pInit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,targetedDim),1));
    tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
end

if includeEarlyInitCohPertTime
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(pInit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,targetedDim),1));
        tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
    end
end

gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
gainRegression = regress(tempData(:,1),[tempData(:,2) ones(size(tempData(:,2)))]);

G = pInit(:,:,:,targetedDim).*gainRegression(1) + gainRegression(2);
G(G<0) = 0;
        
%% Use identified gain over time as gain in feedback model
plotOptsSim.On = false;
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        Gtemp = [G(:,si,ci)'; Gfb(ci)*ones(size(G(:,si,ci)'))];
        input(:,si,ci) = speeds(si)*ones(size(G(:,si,ci)'));
        input(initCoh.neuron_t < 0,si,ci) = 0;
        eta = zeros(size(input(:,si,ci)));
        [eyeSpeed(:,si,ci), C, R] = simulateTightFeedbackFlexibleGains(Gtemp,1,input(:,si,ci)',eta',...
            'assumeSteadyState',false,'plotOpts',plotOptsSim,'latency',80);
    end
end


%% Calculate covariance and correlation over time
C = cov(eyeSpeed);
R = corrcoef(eyeSpeed);

%% Plot results
%% Plot behavioral gain vs activity along targeted dimensions
dimNames = {'Speed','Coherence','Gain','Offset'};
gvthOpt = figure('Name',['Behavioral gain vs activity on optimized targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
plot(tempData(:,2),tempData(:,1),...
    'o','Color',[0.4 0.4 0.4])
hold on

gin = linspace(min(tempData(:,2)),max(tempData(:,2)),20);
plot(gin,gainRegression(1)*gin + gainRegression(2),'r')

axis square
ax = axis;
text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
xlabel('Behavioral gain (unitless)')
ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])

%% Plot responses over time
figure
idx = 1;
for si = 1:length(speeds)
    subplot(1,length(speeds),idx)
    for ci = 1:length(cohs)
        plot(initCoh.neuron_t,input(:,si,ci),'k--')
        hold on
        plot(initCoh.neuron_t,eyeSpeed(:,si,ci),'Color',[speedColors(si,:) cohs(ci)/100])
    end
    idx = idx+1;
end
