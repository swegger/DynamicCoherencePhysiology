function analyzeFitSimpleModelToNeurons(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs

Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
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

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
dirs = Parser.Results.dir;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
resultsFile = Parser.Results.resultsFile;
testInitGain = Parser.Results.testInitGain;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
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
                    ind = find(dirs == initCoh.preferredDirectionRelative(uniti));
                    Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
                    
                    ind = find(dirs == dynCoh.preferredDirectionRelative(uniti));
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
    [Cohs, Spds] = meshgrid(cohsFEF,speedsFEF);
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

%% Find files and get data, fit model parameters, ridge regression value

%% Find model result files
switch subject
    case 'ar'
        files = dir('/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/fitSimpleFEFModelResults/*20240122*.mat');
end

for filei = 1:length(files)
    temp = load([files(filei).folder '/' files(filei).name]);
    if length(temp.modelFEF.tau) == 1
        temp.modelFEF.tau = temp.modelFEF.tau*ones(size(temp.modelFEF.baseLine));
    end
    Rtemp = squeeze(temp.RtoFit_z(:,:,2,:));
    for neuroni = 1:size(Rtemp,3)
        tempHat = SimpleFEFmodel(temp.modelFEF.W(:,neuroni)',...
            temp.modelFEF.baseLine(neuroni),temp.modelFEF.R0(neuroni),temp.modelFEF.tau(neuroni)/temp.modelFEF.dt,temp.inputs_z);
        RhatTemp(:,:,neuroni) = tempHat(:,:,2);
    end
    sse(filei) = sum((Rtemp(:) - RhatTemp(:)).^2);
    sseNeuron(:,filei) = sum( (Rtemp - RhatTemp).^2 ,[1,2]);
    lambdaRidge(filei) = temp.lambdaRidge;
    
end

%% Get data from best fitting model over ridge values
[~,minInd] = min(sse);
temp = load([files(minInd).folder '/' files(minInd).name]);
R = temp.RtoFit_z;
modelFEF = temp.modelFEF;
tau = temp.tau;
dt = temp.dt;
inputs = temp.inputs_z;
for neuroni = 1:size(R,4)
    Rhat(:,:,:,neuroni) = SimpleFEFmodel(modelFEF.W(:,neuroni)',...
        modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),modelFEF.tau(neuroni)/dt,inputs);
end
t = temp.mtResults.mt{1}.neuron_t(temp.iMT);
spref = temp.spref;
spref = spref(~isnan(temp.interpolatedR(1,1,1,:)));
[~,prefSort] = sort(spref);

%%

for si = 1:length(speedsFEF)
    rvalAll = corrcoef([squeeze(R(:,si,2,:)) squeeze(Rhat(:,si,2,:))]);
    rvals(:,si) = diag(rvalAll(size(R,4)+1:end,1:size(R,4)));
end

%% Analysis
inputSigma = squeeze(std(temp.interpolatedR,[],[1,2,3]));
inputSigma = inputSigma(~isnan(temp.interpolatedR(1,1,1,:)));

W = modelFEF.W;
baseLine = modelFEF.baseLine;
if compWithBetaSubspace
    cellID2 = squeeze(temp.cellID(:,1,1:2));
    cellID3 = squeeze(cellID(:,1,1:2));
    LIA = ismember(cellID2,cellID3,'rows');
else
    LIA = true(size(cellID2,1));
end
W = W(:,LIA);
baseLine = baseLine(LIA);
Wz = W.*repmat(inputSigma,[1,size(W,2)]);

% Normalize by maximium value of W for each FEF neuron
wmax = max(abs(Wz),[],2);
Weff = Wz./repmat(wmax,[1,size(W,2)]);

% SVD
[U,S,V] = svd(W);

% Approach from the important dimensions in Rhat
Y = nan(size(R,1)*size(R,2)*size(R,3),sum(LIA));
Yhat = nan(size(Rhat,1)*size(Rhat,2)*size(Rhat,3),sum(LIA));
idx = 1;
for si = 1:size(Rhat,2)
    for ci = 1:size(Rhat,3)
        Y(idx:idx+size(R,1)-1,:) = R(:,si,ci,LIA);
        Yhat(idx:idx+size(Rhat,1)-1,:) = Rhat(:,si,ci,LIA);
        idx = idx+size(Rhat,1);
    end
end
[Uy,Sy,Vy] = svd(Y,0);
Wy = W*Vy;
[Uhat,Shat,Vhat] = svd(Yhat,0);
What = W*Vhat;                      % Now ordering weight matrix by the important dimensions in Rhat (e.g. reduced rank regression approach)

[~,rvalsort] = sort(rvals(LIA,3));
%% Plot results

%% Best ridge regression for predicting held-out data
figure
[~,sortInd] = sort(lambdaRidge);
subplot(1,2,1)
loglog(lambdaRidge(sortInd),sse(sortInd),'ko-')
hold on
plotVertical(lambdaRidge(sse == min(sse)));
xlabel('Ridge regression coefficient')
ylabel('Sum of squared errors (spikes/sec)^2')

subplot(1,2,2)
loglog(lambdaRidge(sortInd),sseNeuron(:,sortInd),'-','Color',[0 0 0 0.1])
[~,sseMins] = min(sseNeuron,[],2);
hold on
plotVertical(lambdaRidge(sseMins));
xlabel('Ridge regression coefficient')
ylabel('Sum of squared errors (spikes/sec)^2')


%% Best and worst fits (on average across held-out data)
mrvals = mean(rvals,2);
[~,mrsort] = sort(mrvals);

figure
for ni = 1:5
    subplot(2,5,ni)
    plot(t,R(:,:,2,mrsort(ni))*1000)
    hold on
    plot(t,Rhat(:,:,2,mrsort(ni))*1000,'k')
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    title(['Rank = ' num2str(ni)])
end

for ni = 1:5
    subplot(2,5,ni+5)
    plot(t,R(:,:,2,mrsort(end-ni+1))*1000)
    hold on
    plot(t,Rhat(:,:,2,mrsort(end-ni+1))*1000,'k')
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    title(['Rank = ' num2str(length(mrsort)-ni+1)])
end

%% Prediction of held-out data
figure
subplot(1,3,1)
title('Example neuron')
neuroni = find(ismember(cellID2, [77,140],'rows'));
for si = 1:length(speedsFEF)
    plot(t,R(:,si,2,neuroni)*1000,'Color',speedColors(si,:))
    hold on
    plot(t,Rhat(:,si,2,neuroni)*1000,'--','Color',speedColors(si,:))
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

subplot(1,3,2)
title('Mean across neurons')
for si = 1:length(speedsFEF)
    plot(t,nanmean(R(:,si,2,:),4)*1000,'Color',speedColors(si,:))
    hold on
    plot(t,nanmean(Rhat(:,si,2,:),4)*1000,'--','Color',speedColors(si,:))
end
xlabel('Time from motion onset (ms)')
ylabel('Mean spikes/s (across neurons)')

subplot(1,3,3)
title('R^2 by neuron')
for si = 1:length(speedsFEF)
    plot(rvals(:,si).^2,'o-','Color',speedColors(si,:))
    hold on
end
xlim([1 256])
xlabel('Neuron number')
ylabel('R^2')

%% Model parameters
figure
ah = subplot(1,2,1);
imagesc(W(prefSort,rvalsort)')
ah.YDir = 'normal';
xlabel('MT neuron')
ylabel('FEF neuron')

subplot(1,2,2)
plot(baseLine(rvalsort)*1000,1:length(baseLine),'o')
axis tight
ylabel('FEF  neuron')
xlabel('Baseline rate')

%% Analysis of weight matrix
figure
subplot(2,2,1)
imagesc(W(prefSort,rvalsort)')
xlabel('MT neuron')
ylabel('FEF neuron')

subplot(2,2,2)
plot(cumsum(diag(S.^2))./sum(diag(S.^2)),'o')
axis tight
xlabel('Dimensions')
ylabel('Cumulative variance explained')

subplot(2,2,3)
rankN = 1;
for ri = 1:rankN
    semilogx(spref(prefSort),U(prefSort,ri),'o')
    hold on
end
xlim([1 300])
xlabel('Speed preference (deg/s)')
ylabel('MT unit weight')

subplot(2,2,4)
for ri = 1:rankN
    plot(V(rvalsort,ri),'o')
    hold on
end
xlabel('FEF neuron #')
ylabel('FEF unit weight')

%% Analysis of weight matrix in space that best describes predicted responses
figure
subplot(1,3,1)
plot(cumsum(diag(Shat.^2))./sum(diag(Shat.^2)),'o')
axis tight
xlabel('Dimensions')
ylabel('Cumulative variance explained')

subplot(1,3,2)
rankN = 1;
for ri = 1:rankN
    semilogx(spref(prefSort),What(prefSort,ri),'o')
    hold on
end
xlim([1 300])
xlabel('Speed preference (deg/s)')
ylabel('MT unit weight')

subplot(1,3,3)
for ri = 1:rankN
    plot(Vhat(rvalsort,ri),'o')
    hold on
end
xlabel('FEF neuron #')
ylabel('FEF unit weight')

figure
ci = 2;
for ri = 1:rankN
    subplot(rankN,4,1+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,-(squeeze(R(:,si,ci,LIA))*Vhat(:,ri) - baseLine*Vhat(:,ri)),'Color',speedColors(si,:))
        hold on
        plot(t,-(squeeze(Rhat(:,si,ci,LIA))*Vhat(:,ri) - baseLine*Vhat(:,ri)),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Decay dynamics along dim ' num2str(ri)])
    
    subplot(rankN,4,2+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,What(:,ri)'*inputs(:,:,si,ci),'Color',speedColors(si,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['MT inputs along dim ' num2str(ri)])
    
    subplot(rankN,4,3+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,-(squeeze(R(:,si,ci,LIA))*Vhat(:,ri) - baseLine*Vhat(:,ri))' + What(:,ri)'*inputs(:,:,si,ci),'Color',speedColors(si,:))
        hold on
        plot(t,-(squeeze(Rhat(:,si,ci,LIA))*Vhat(:,ri) - baseLine*Vhat(:,ri))' + What(:,ri)'*inputs(:,:,si,ci),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['dr/dt along dim ' num2str(ri)])
    
    subplot(rankN,4,4+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,squeeze(R(:,si,ci,LIA))*Vhat(:,ri),'Color',speedColors(si,:))
        hold on
        plot(t,squeeze(Rhat(:,si,ci,LIA))*Vhat(:,ri),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Response along dim ' num2str(ri)])
    
end

%% Analysis of weight matrix in space that best describes speed/coherence/gain
figure
Wsubspace = W*BinitOrth(:,1:4)*BinitOrth(:,1:4)';
baseLineSubspace = baseLine*BinitOrth(:,1:4)*BinitOrth(:,1:4)';
for dimi = 1:size(BinitTi,1)
    subplot(size(BinitTi,1)+1,2,1+(dimi-1)*2)
    WsubspaceR = W*BinitOrth(:,dimi)*BinitOrth(:,dimi)';
    baseLineSubspaceR = baseLine*BinitOrth(:,dimi)*BinitOrth(:,dimi)';
    imagesc(WsubspaceR(prefSort,rvalsort)')
    set(gca,'YDir','normal')
    axis tight
    xlabel('MT neuron (sorted by speed preference)')
    ylabel('FEF neuron (sorted by rvalue of fit)')
    title(['Dimension #' num2str(dimi)])
    
    subplot(size(BinitTi,1)+1,2,2+(dimi-1)*2)
    plot(baseLineSubspaceR(rvalsort),1:size(W,2),'o')
    axis tight
    ylabel('Neuron #')
    xlabel('Baseline along targeted dimensions')
end
dimi = dimi+1;
subplot(size(BinitTi,1)+1,2,1+(dimi-1)*2)
imagesc(Wsubspace(prefSort,rvalsort)')
set(gca,'YDir','normal')
axis tight
xlabel('Dimensions')
ylabel('Cumulative variance explained')

subplot(size(BinitTi,1)+1,2,2+(dimi-1)*2)
plot(baseLineSubspace(rvalsort),1:size(W,2),'o')
axis tight
ylabel('Neuron #')
xlabel('Baseline along targeted dimensions')

figure
for dimi = 1:size(BinitTi,1)
    subplot(size(BinitTi,1),2,1+(dimi-1)*2)
    WsubspaceR = W*BinitOrth(:,dimi);
    baseLineSubspaceR = baseLine*BinitOrth(:,dimi)*BinitOrth(:,dimi)';
    semilogx(spref,WsubspaceR(prefSort),'o')
    axis tight
    xlim([0.1 128])
    xlabel('Speed preference (deg/s)')
    ylabel(['Weight along dimension #' num2str(dimi)])
    title(['Dimension #' num2str(dimi)])
    
    subplot(size(BinitTi,1),2,2+(dimi-1)*2)
    plot(rvalsort,BinitTi(dimi,rvalsort),'o')
    axis tight
    xlabel('FEF Neuron # (sorted by rvalue of model fit)')
    ylabel(['Weight along dimension #' num2str(dimi)])
end



figure
ci = 2;
Wsub = W*BinitOrth;
for ri = 1:size(BinitTi,1)
    subplot(size(BinitTi,1),4,1+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,-(squeeze(R(:,si,ci,LIA))*BinitOrth(:,ri) - baseLine*BinitOrth(:,ri)),'Color',speedColors(si,:))
        hold on
        plot(t,-(squeeze(Rhat(:,si,ci,LIA))*BinitOrth(:,ri) - baseLine*BinitOrth(:,ri)),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Decay dynamics along dim ' num2str(ri)])
    
    subplot(size(BinitTi,1),4,2+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,Wsub(:,ri)'*inputs(:,:,si,ci),'Color',speedColors(si,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['MT inputs along dim ' num2str(ri)])
    
    subplot(size(BinitTi,1),4,3+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,-(squeeze(R(:,si,ci,LIA))*BinitOrth(:,ri) - baseLine*BinitOrth(:,ri))' + Wsub(:,ri)'*inputs(:,:,si,ci),'Color',speedColors(si,:))
        hold on
        plot(t,-(squeeze(Rhat(:,si,ci,LIA))*BinitOrth(:,ri) - baseLine*BinitOrth(:,ri))' + Wsub(:,ri)'*inputs(:,:,si,ci),'--','Color',speedColors(si,:))   
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['dr/dt along dim ' num2str(ri)])
    
    subplot(size(BinitTi,1),4,4+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,squeeze(R(:,si,ci,LIA))*BinitOrth(:,ri),'Color',speedColors(si,:))
        hold on
        plot(t,squeeze(Rhat(:,si,ci,LIA))*BinitOrth(:,ri),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Response along dim ' num2str(ri)])
    
end

%% Functions

%% Simple model
function R = SimpleFEFmodel(W,baseLine,R0,tau,inputs)
    sz = size(inputs);
    R = nan(sz(2:4));
    R(1,:,:) = repmat(R0,[1,sz(3),sz(4)]);
    for ti = 2:sz(2)
        for si = 1:sz(3)
            for ci = 1:sz(4)
                dR(ti,si,ci) = -(R(ti-1,si,ci)-baseLine) + W*inputs(:,ti,si,ci);
                R(ti,si,ci) = R(ti-1,si,ci) + dR(ti,si,ci)/tau;
            end
        end
    end