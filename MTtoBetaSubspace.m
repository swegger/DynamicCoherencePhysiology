function MTtoBetaSubspace(varargin)
%% 
%
%
%
%%

%% Defaults
speedPrefOpts_default.tWin = [40,120];
speedPrefOpts_default.P0 = [16,1];
speedPrefOpts_default.ub = [254,128];
speedPrefOpts_default.lb = [0, 0];
speedPrefOpts_default.c = NaN;
speedPrefOpts_default.s = NaN;
speedPrefOpts_default.d = 0;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speedsMT',[2,4,8,16,32])
addParameter(Parser,'cohsMT',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
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
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'P0',NaN)
addParameter(Parser,'lb',NaN)
addParameter(Parser,'ub',NaN)

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
objectFile = Parser.Results.objectFile;
speedsMT = Parser.Results.speedsMT;
cohsMT = Parser.Results.cohsMT;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
dir = Parser.Results.dir;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
resultsFile = Parser.Results.resultsFile;
testInitGain = Parser.Results.testInitGain;
includeEarlyInitCohPertTime = Parser.Results.includeEarlyInitCohPertTime;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
sequences = Parser.Results.sequences;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
sprefFromFit = Parser.Results.sprefFromFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
checkMTFit = Parser.Results.checkMTFit;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
P0 = Parser.Results.P0;
lb = Parser.Results.lb;
ub = Parser.Results.ub;


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

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
[MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
    'speedsMT',speedsMT,'cohsMT',cohsMT,'directionsMT',directionsMT,...
    'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
    'speedPrefOpts',speedPrefOpts);

% MT = [];
% spref = [];
% swidth = [];
% dirMT = [];
% if opponentMT
%     dN = 1;
% else
%     dN = length(directionsMT);
% end
% for di = 1:dN
%     for filei = 1:length(mtResults.mt)
%         disp(['File ' num2str(filei) ' of ' num2str(length(mtResults.mt))])
%     
%         MTtemp = nan(length(mtResults.mt{filei}.neuron_t),length(speedsMT),length(cohsMT));
%         
%         
%         for si = 1:length(speedsMT)
%             for ci = 1:length(cohsMT)
%                 if opponentMT
%                     [~,condLogical] = trialSort(mtResults.mt{filei},directionsMT(1),speedsMT(si),NaN,cohsMT(ci));
%                     [~,condLogicalNull] = trialSort(mtResults.mt{filei},directionsMT(2),speedsMT(si),NaN,cohsMT(ci));
%                     MTtemp(:,si,ci) = mean(mtResults.mt{filei}.r(:,condLogical),2) - mean(mtResults.mt{filei}.r(:,condLogicalNull),2);
%                 else
%                     [~,condLogical] = trialSort(mtResults.mt{filei},directionsMT(di),speedsMT(si),NaN,cohsMT(ci));
%                     MTtemp(:,si,ci) = mean(mtResults.mt{filei}.r(:,condLogical),2);
%                 end
%             end
%         end
%         
%         if sprefFromFit
%             [mu,sig,~,normR,~] = fitSpeedTuning(mtResults. mt{filei},'P0',speedPrefOpts.P0,...
%                 'ub',speedPrefOpts.ub,'lb',speedPrefOpts.lb,...
%                 'c',speedPrefOpts.c,'s',speedPrefOpts.s,'d',speedPrefOpts.d,...
%                 'tWin',speedPrefOpts.tWin);
%             spref = [spref mu];
%             swidth = [swidth sig];
%             
%             if checkMTFit
%                 s = linspace(min(speedsMT)-0.1*min(speedsMT),max(speedsMT)*1.1,20);
%                 h = figure;
%                 semilogx(mtResults.mt{filei}.speeds,normR,'o')
%                 hold on
%                 semilogx(s,speedTuning(mtResults.mt{filei},s,mu,sig))
%                 ax = axis;
%                 text(ax(1),0.95*ax(4),['\mu = ' num2str(mu) ', \sig = ' num2str(sig)])
%                 input('Press enter to continue ')
%                 close(h);
%             end
%         else
%             spref = [spref speedsMT(nansum(MTtemp(mtResults.mt{filei}.neuron_t >= 50 & mtResults.mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
%                 max(nansum(MTtemp(mtResults.mt{filei}.neuron_t >= 50 & mtResults.mt{filei}.neuron_t <= 150,:,cohs == 100))))];
%         end
%         dirMT = [dirMT, directionsMT(di)];
%         MT = cat(4,MT,MTtemp);
%     end
% end

mtNeuron_t = mtResults.mt{filei}.neuron_t;


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speedsMT,cohsMT,speeds,cohs);

inputs = permute(interpolatedR,[4,1,2,3]);
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);

%% Iterate across neurons and fit data

fef_t = initCoh.neuron_t(initCoh.neuron_t>=-100 & initCoh.neuron_t<=1350);
% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,mtNeuron_t);
RtoFit = pInit(iFEF,:,:,:);
RtoFit_z = (RtoFit-mean(RtoFit,[1,2,3]))./std(RtoFit,[],[1,2,3]);
inputs = inputs(:,iMT,:,:);
dt = mtNeuron_t(2)-mtNeuron_t(1);
t = mtNeuron_t(iMT);
inputs_z = (inputs-mean(inputs(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(inputs(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);

%% Regression models

ridgeLambda = logspace(-1,12,10);

% Build matrices
Xfit = nan(size(inputs,2)*size(inputs,3)*(size(inputs,4)-1),size(inputs,1)+1);
Yfit = nan(size(RtoFit,1)*size(RtoFit,2)*(size(RtoFit,3)-1),size(RtoFit,4));

idx = [1, 1];
for si = 1:length(speeds)
    for ci = [1 length(cohs)]
        Xfit(idx(1):idx(1)+size(inputs,2)-1,:) = [permute(inputs_z(:,:,si,ci),[2,1]) ones(size(inputs,2),1)];
        Yfit(idx(2):idx(2)+size(RtoFit,1)-1,:) = squeeze(RtoFit_z(:,si,ci,:));
        idx(1) = idx(1) + size(inputs,2);
        idx(2) = idx(2) + size(RtoFit,1);
    end
end

Xtest = nan(size(inputs,2)*size(inputs,3),size(inputs,1)+1);
Ytest = nan(size(RtoFit,1)*size(RtoFit,2),size(RtoFit,4));
idx = [1, 1];
for si = 1:length(speeds)
    Xtest(idx(1):idx(1)+size(inputs,2)-1,:) = [permute(inputs_z(:,:,si,2),[2,1]) ones(size(inputs,2),1)];
    Ytest(idx(2):idx(2)+size(RtoFit,1)-1,:) = squeeze(RtoFit_z(:,si,2,:));
    idx(1) = idx(1) + size(inputs,2);
    idx(2) = idx(2) + size(RtoFit,1);
end

% Regress
trainBeta = nan([size(Xfit,2),size(Yfit,2),length(ridgeLambda)]);
for ri = 1:length(ridgeLambda)
    trainBeta(:,:,ri) = inv(Xfit'*Xfit + ridgeLambda(ri)*eye(size(Xfit'*Xfit)))*Xfit'*Yfit;
    trainYfit = Xfit*trainBeta(:,:,ri);
    testYfit = Xtest*trainBeta(:,:,ri);
    
    fitQ(ri,1) = sum( (Yfit(:) - trainYfit(:)).^2 );
    fitQ(ri,2) = sum( (Ytest(:) - testYfit(:)).^2 );
end
Beta = trainBeta(:,:,fitQ(:,2)==min(fitQ(:,2)));

% Predicted responses
YhatTrain = Xfit*Beta;
YhatTest = Xtest*Beta;

[Uhat,Shat,Vhat] = svd(YhatTrain,0);
BetaBar = Beta*Vhat;

Rhat = nan(size(RtoFit));
idx = 1;
idxTest = 1;
for si = 1:length(speeds)
    for ci = [1,length(cohs)]
        Rhat(:,si,ci,:) = YhatTrain(idx:idx+size(RtoFit,1)-1,:);
        idx = idx+size(RtoFit,1);
    end
    
    ci = 2;
    Rhat(:,si,ci,:) = YhatTest(idxTest:idxTest+size(RtoFit,1)-1,:);
    idxTest = idxTest+size(RtoFit,1);
end

% Find correlation with held out data and do Fisher transformation
rvals = corrcoef([Ytest YhatTest]);
RhatCC = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZ = 0.5*log( (1+RhatCC)./(1-RhatCC) );

spref2 = spref(~isnan(interpolatedR(1,1,1,:)));

    % Fit quality
    
    t = fef_t(iFEF);
 
    figure
    for neuroni = 1:size(RtoFit,4)
        for ci = 1:length(cohs)
            subplot(size(RtoFit,4),length(cohs),ci + (neuroni-1)*length(cohs))
            for si = 1:length(speeds)
                plot(t,RtoFit_z(:,si,ci,neuroni),'Color',speedColors(si,:))
                hold on
                plot(t,Rhat(:,si,ci,neuroni),'--','Color',speedColors(si,:))
            end
            xlabel('Time from motion onset (ms)')
            ylabel(['Response along dim ' num2str(neuroni)])
        end
    end
    
    % Weighted MT inputs
    [~,sortInd] = sort(spref2);
    cax = [Inf,-Inf];
    figure
    for axi = 1:3
        for si = 1:length(speeds)
            subplot(3,length(speeds),si+(axi-1)*length(speeds))
            temp = inputs_z(:,:,si,3).*Beta(1:end-1,axi); %log2(spref2');%
            imagesc(t,1:length(spref2),temp(sortInd,:))
            caxTemp = caxis;
            cax(1) = min([cax(1) caxTemp(1)]);
            cax(2) = max([cax(2) caxTemp(2)]);
        end
    end
    for subploti = 1:3*length(speeds)
        subplot(3,length(speeds),subploti)
        caxis(cax)
        xlabel('Time from motion onset (ms)')
        ylabel('MT neuron # (sorted by speed preference)')
    end
    
%% Fit dynamic model

lambdaRidge = ridgeLambda(fitQ(:,2)==min(fitQ(:,2)));

if any(isnan(P0))
    P0 = [20/dt,0,randn(1,size(inputs,1))/1000];
end
if any(isnan(lb))
    lb = [1, -Inf, -Inf*ones(1,size(inputs,1))];
end
if any(isnan(ub))
    ub = [Inf, 100, Inf*ones(1,size(inputs,1))];
end

OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10);
for neuroni = 1:size(RtoFit,4)
    tic
    disp(['Dim ' num2str(neuroni) ' of ' num2str(size(RtoFit,4))])
    Rtemp = RtoFit_z(:,:,[1, length(cohs)],neuroni);
    R0 = nanmean(Rtemp(1,:,:),[2,3]);
%     if neuroni == 3
%         ztemp = cat(1,permute(RtoFit(:,:,:,1),[4,1,2,3]),permute(RtoFit(:,:,:,2),[4,1,2,3]));
%         ztemp = (ztemp-mean(ztemp,[2,3,4]))./std(ztemp,[],[2,3,4]);
%         minimizer = @(P)minimizant(P,R0,ztemp(:,:,:,[1,length(cohs)]),Rtemp,0);
%         [Ptemp, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,[P0(1:2) 0 0],[],[],[],[],[lb(1:2) -Inf -Inf],[ub(1:2) Inf Inf],[],OPTIONS);
%         P = [Ptemp(1:end-2) zeros(1,size(inputs,1))];
%         A = Ptemp(end-1:end);
%         grad = zeros(size(P0,2),1);
%         hessian = zeros(size(P0,2),size(P0,2));
%     else
        minimizer = @(P)minimizant(P,R0,inputs_z(:,:,:,[1,length(cohs)]),Rtemp,lambdaRidge);
        [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
        A = [0; 0];
%     end
    model.tau(neuroni) = P(1)*dt;
    model.R0(neuroni) = R0;
    model.baseLine(neuroni) = P(2);
    model.W(:,neuroni) = P(3:end);
    model.A(:,neuroni) = A;
    model.fval(:,neuroni) = fval;
    model.exitflag(:,neuroni) = exitflag;
    model.output(:,neuroni) = output;
    model.lambda(:,neuroni) = lambda;
    model.grad(:,neuroni) = grad;
    model.hessian(:,:,neuroni) = hessian;
    toc
end

model.dt = dt;


%% Plotting

%% Fit quality
    t = fef_t(iFEF);
 
    figure
    for neuroni = 1:size(RtoFit,4)
        
        transformBeta = regress(model.W(:,neuroni),[log2(spref2') ones(size(model.W(:,neuroni)))]);
        tempFEF = RtoFit_z(:,:,:,neuroni);
%         if neuroni == 3
%             ztemp = cat(1,permute(RtoFit(:,:,:,1),[4,1,2,3]),permute(RtoFit(:,:,:,2),[4,1,2,3]));
%             ztemp = (ztemp-mean(ztemp,[2,3,4]))./std(ztemp,[],[2,3,4]);
%             test = SimpleFEFmodel(model.A(:,neuroni)',model.baseLine(neuroni),model.R0(neuroni),model.tau(neuroni)/model.dt,...
%                 ztemp);
%         else

%             test = SimpleFEFmodel(([log2(spref2') ones(size(model.W(:,axi)))]*transformBeta)',model.baseLine(neuroni),model.R0(neuroni),model.tau(neuroni)/model.dt,inputs_z);
            test = SimpleFEFmodel(model.W(:,neuroni)',model.baseLine(neuroni),model.R0(neuroni),model.tau(neuroni)/model.dt,inputs_z);
%         end
        for ci = 1:length(cohs)
            subplot(size(RtoFit,4),length(cohs),ci + (neuroni-1)*length(cohs))
            for si = 1:length(speeds)
                plot(t,tempFEF(:,si,ci),'Color',speedColors(si,:))
                hold on
                plot(t,test(:,si,ci),'--','Color',speedColors(si,:))
            end
            xlabel('Time from motion onset (ms)')
            ylabel(['Response along dim ' num2str(neuroni)])
        end
    end

%% Weight matrix analysis
figure
for axi = 1:size(RtoFit,4)
    subplot(1,size(RtoFit,4),axi)
    plot(log2(spref2),model.W(:,axi),'o')
end

%% Weighted inputs
figure
[~,sortInd] = sort(spref2);
[~,sortIndW] = sort(model.W(:,1));
    cax = [Inf,-Inf];
    for axi = 1:3
        transformBeta = regress(model.W(:,axi),[log2(spref2') ones(size(model.W(:,axi)))]);
        for si = 1:length(speeds)
            subplot(3,length(speeds),si+(axi-1)*length(speeds))
%             temp = inputs_z(:,:,si,3).*([log2(spref2') ones(size(model.W(:,axi)))]*transformBeta);
            temp = inputs_z(:,:,si,3).*model.W(:,axi);% - inputs_z(:,:,si,3).*([log2(spref2') ones(size(model.W(:,axi)))]*transformBeta);
            imagesc(t,1:length(spref2),temp(sortInd,:))
            caxTemp = caxis;
            cax(1) = min([cax(1) caxTemp(1)]);
            cax(2) = max([cax(2) caxTemp(2)]);
        end
    end
    for subploti = 1:3*length(speeds)
        subplot(3,length(speeds),subploti)
        caxis(cax)
        xlabel('Time from motion onset (ms)')
        ylabel('MT neuron # (sorted by speed preference)')
    end

%% Plot targeted dimensionality reduction
figure('Name','Input related activity','Position',[63 169 1606 1079])

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
            pInit(initTaccept,speedi,cohi,4),pInit(initTaccept,speedi,cohi,3),...
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
        pDynCross(dynTaccept,seqi,4),pDynCross(dynTaccept,seqi,3),...
        '-','Color',colors(seqi,:))
    hold on
end
grid on
xlabel('Speed related activity')
ylabel('Gain related activity')

subplot(1,3,3)
for seqi = 1:5
    plot(pDyn(dynTaccept,seqi,3),pDyn(dynTaccept,seqi,2),'-','Color',colors(seqi,:))
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


%% Functions

%% Function to fit
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
    
 %% Objective
function out = minimizant(P,R0,inputs,R,lambda)
    tau = P(1);
    baseLine = P(2);
    W = P(3:end);
    Rest = SimpleFEFmodel(W,baseLine,R0,tau,inputs);
    out = sum((R(:) - Rest(:)).^2) + lambda*sum(W.^2);