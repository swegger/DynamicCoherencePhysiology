function Beta = MTtoFEFregression(dcp,varargin)
%% MTtoFEFregression
%
%
%%

%% Defaults
plotOpts_default.On = false;

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
compToBehavioralGain_default.On = false;
neuronTypingComparison_default.On = false;

theoretical_default.weightTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';

simulateMT_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'dynCohCollate',false)
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'trainCondition',trainCondition_default)
addParameter(Parser,'tWin',[0 900])
addParameter(Parser,'rankN',80)
addParameter(Parser,'ridgeLambda',logspace(-1,12,10))
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'simulateMT',simulateMT_default)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'compToBehavioralGain',compToBehavioralGain_default)
addParameter(Parser,'neuronTypingComparison',neuronTypingComparison_default)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'plotOpts',plotOpts_default)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
objectFile = Parser.Results.objectFile;
dynCohCollate = Parser.Results.dynCohCollate;
chanMap = Parser.Results.chanMap;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
trainCondition = Parser.Results.trainCondition;
directions = Parser.Results.directions;
tWin = Parser.Results.tWin;
rankN = Parser.Results.rankN;
ridgeLambda = Parser.Results.ridgeLambda;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
simulateMT = Parser.Results.simulateMT;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
compToBehavioralGain = Parser.Results.compToBehavioralGain;
neuronTypingComparison = Parser.Results.neuronTypingComparison;
theoretical = Parser.Results.theoretical;
rateCutoff = Parser.Results.rateCutoff;
checkUnitType = Parser.Results.checkUnitType;
plotOpts = Parser.Results.plotOpts;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

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

[RinitTemp, RdynTemp, cellIDTemp, passCutoffTemp, locationsTemp] = collateFiringRates(dcp,...
    'sourceDirectory',sourceDirectory,'directions',directions,'chanMap',chanMap,...
    'rateCutoff',rateCutoff,'checkUnitType',checkUnitType,...
    'initCohCollate',true,'dynCohCollate',dynCohCollate);

% Add to data matrices
Rinit(:,:,:,indx:indx+size(RinitTemp,4)-1) = RinitTemp;
Rdyn(:,:,indx:indx+size(RdynTemp,3)-1) = RdynTemp;
cellID(indx:indx+size(cellIDTemp,1)-1,:,1:3) = cellIDTemp;
passCutoff(indx:indx+size(passCutoffTemp,1)-1,:) = passCutoffTemp;
locations(indx:indx+size(locationsTemp,1)-1,:) = locationsTemp;

indx = indx+size(locationsTemp,1);

temp = load([sourceDirectory '/' dcp{1}.datapath(end-8:end-1) 'obj/initCoh' ...
    dcp{1}.datapath(end-8:end)]);
initCoh = temp.initCoh;

if dynCohCollate
    filei = 1;
    successfulDynCohLoad = false;
    while ~successfulDynCohLoad
        tempFile = [sourceDirectory '/' dcpDynCoh.dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
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
% indx = 1;
% for filei = 1:length(dcp)
%     disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
%     
%     % Add probe info
%     dcp{filei} = addProbeInfo(dcp{filei});
%     
%     % InitCoh data
%     load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
%         dcp{filei}.datapath(end-8:end)])
%     
%     
%     if ~isempty(initCoh.R)
%         
%         passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;
%         
%         % Get data for each neuron
%         for uniti = 1:length(initCoh.preferredDirectionRelative)
%             ind = find(directions == initCoh.preferredDirectionRelative(uniti));
%             Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
%             
%             
%             for j = 1:length(initCoh.unitIndex)
%                 cellID(indx,j,1) = filei;
%                 cellID(indx,j,2) = initCoh.unitIndex(uniti);
%                 cellID(indx,j,3) = initCoh.unitIndex(j);
%             end
%                         
%             indx = indx+1;
%         end
%     end
% end
% 
% 
% Rinit = Rinit(:,:,:,1:indx-1);
% passCutoff = logical(passCutoff(1:indx-1));
% cellID = cellID(1:indx-1,:,:);

%taccept = initCoh.neuron_t >= tWin(1) & initCoh.neuron_t <= tWin(2);

%Rinit = Rinit(taccept,:,:,:);
fef_t = initCoh.neuron_t;%(taccept);

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
if dynCohCollate
    m2 = squeeze(max(Rdyn,[],[1,2]))*1000;
else
    m2 = zeros(size(squeeze(max(Rdyn,[],[1,2]))*1000));
end
Rinit = Rinit(:,:,:,m<=150 & m2<=150);
Rdyn = Rdyn(:,:,m<=150 & m2<=150);
locations = locations(m<=150 & m2<=150,:);
cellID = cellID(m<=150 & m2<150,:,:);

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
% Get mean of each unit
if simulateMT.On
    mtNeuron_t = -500:900;
    speeds = speedsFEF;
    [MT, spref, cohs] = simulateMTdata(mtResults,mtNeuron_t,simulateMT.modelN,...
        'speedsMT',speeds,'removeBaseline',simulateMT.removeBaseline,'gaussianApprox',simulateMT.gaussianApprox,...
        'accept',simulateMT.accept);
    swidth = 1.2*ones(size(spref));
else
    [MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
        'speedsMT',speeds,'cohsMT',cohs,'directionsMT',directionsMT,...
        'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
        'speedPrefOpts',speedPrefOpts);
    
    mtNeuron_t = mtResults.mt{1}.neuron_t;
end

% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);

inputs = permute(interpolatedR,[4,1,2,3]);
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);

% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,mtNeuron_t);
R = Rinit(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);
t = mtNeuron_t(iMT);
inputs_z = (inputs-mean(inputs(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(inputs(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);
R_z = (R-mean(R,[1,2,3]))./std(R,[],[1,2,3]);
taccept = t >= tWin(1) & t <= tWin(2);
inputsToFit = inputs_z(:,taccept,:,:);
RtoFit = R(taccept,:,:,:);
RtoFit_z = (RtoFit-mean(RtoFit,[1,2,3]))./std(RtoFit,[],[1,2,3]);

%% Regression models

% Build matrices
Xfit = nan(size(inputsToFit,2)*sum(trainCondition(:)),size(inputsToFit,1));
Yfit = nan(size(RtoFit_z,1)*sum(trainCondition(:)),size(RtoFit_z,4));
Xall = nan(size(inputs,2)*numel(trainCondition),size(inputs,1));
conditions = nan(size(RtoFit_z,1)*sum(trainCondition(:)),3);
conditionsAll = nan(size(R_z,1)*sum(trainCondition(:)),2);

idx = [1, 1];
idxAll = [1, 1];
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        if trainCondition(si,ci)
            Xfit(idx(1):idx(1)+size(inputsToFit,2)-1,:) = permute(inputsToFit(:,:,si,ci),[2,1]);
            Yfit(idx(2):idx(2)+size(RtoFit_z,1)-1,:) = squeeze(RtoFit_z(:,si,ci,:));
            conditions(idx(2):idx(2)+size(RtoFit_z,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit_z,1),1]);
            idx(1) = idx(1) + size(inputsToFit,2);
            idx(2) = idx(2) + size(RtoFit_z,1);
        end
        Xall(idxAll(1):idxAll(1)+size(inputs,2)-1,:) = permute(inputs_z(:,:,si,ci),[2,1]);
        conditionsAll(idxAll(2):idxAll(2)+size(R_z,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci)],[size(R_z,1),1]);
        idxAll(1) = idxAll(1) + size(inputs,2);
        idxAll(2) = idxAll(2) + size(R_z,1);
    end
end

Xtest = nan(size(inputsToFit,2)*sum(~trainCondition(:)),size(inputsToFit,1));
Ytest = nan(size(RtoFit_z,1)*sum(~trainCondition(:)),size(RtoFit_z,4));
conditionsTest = nan(size(RtoFit_z,1)*sum(~trainCondition(:)),3);
idx = [1, 1];
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        if ~trainCondition(si,ci)
            Xtest(idx(1):idx(1)+size(inputsToFit,2)-1,:) = permute(inputsToFit(:,:,si,2),[2,1]);
            Ytest(idx(2):idx(2)+size(RtoFit_z,1)-1,:) = squeeze(RtoFit_z(:,si,2,:));
            conditionsTest(idx(2):idx(2)+size(RtoFit_z,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit_z,1),1]);
            idx(1) = idx(1) + size(inputsToFit,2);
            idx(2) = idx(2) + size(RtoFit_z,1);
        end
    end
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
Yhat = Xall*Beta;

% Reduced rank regression
[Uhat,Shat,Vhat] = svd(YhatTrain,0);
for ri = 1:rankN
    BetaRRR = Beta*Vhat(:,1:ri)*Vhat(:,1:ri)';
    YrrrFit = Xfit*BetaRRR;
    YrrrTest = Xtest*BetaRRR;
    ccTemp = corrcoef(YrrrFit(:),Yfit(:));
    fitRank(ri,1) = ccTemp(1,2);
    ccTemp = corrcoef(YrrrTest(:),Ytest(:));
    fitRank(ri,2) = ccTemp(1,2);
    
%     fitRank(ri,1) = sum( (YrrrFit(:) - Yfit(:)).^2 );
%     fitRank(ri,2) = sum( (YrrrTest(:) - Ytest(:)).^2 );
end
[~,rankOpt] = max(fitRank(:,2));

BetaRRR = Beta*Vhat(:,1:rankOpt)*Vhat(:,1:rankOpt)';
BetaNull = Beta*Vhat(:,end-rankOpt+1:end)*Vhat(:,end-rankOpt+1:end)';
%YrrrTrain = Xfit*BetaRRR;
%YnullTrain = Xfit*BetaNull;
YrrrTest = Xtest*BetaRRR;
YnullTest = Xtest*BetaNull;
Yrrr = Xall*BetaRRR;
Ynull = Xall*BetaNull;

Rhat = nan(size(R_z));
Rnull = nan(size(R_z));
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        temp = ismember(conditionsAll,[speedsFEF(si),cohsFEF(ci)],'rows');
        Rhat(:,si,ci,:) = Yrrr(temp,:);
        Rnull(:,si,ci,:) = Ynull(temp,:);
%         if trainCondition(si,ci)
%             temp = ismember(conditions,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
%             Rhat(:,si,ci,:) = YrrrTrain(temp,:);
%             Rnull(:,si,ci,:) = YnullTrain(temp,:);
%         else
%             temp = ismember(conditionsTest,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
%             Rhat(:,si,ci,:) = YrrrTest(temp,:);
%             Rnull(:,si,ci,:) = YnullTest(temp,:);
%         end
    end
end

% Find correlation with held out data and do Fisher transformation
rvals = corrcoef([Ytest YhatTest]);
RhatCC_full = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZ_full = 0.5*log( (1+RhatCC_full)./(1-RhatCC_full) );
rvals = corrcoef([Ytest YrrrTest]);
RhatCC = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZ = 0.5*log( (1+RhatCC)./(1-RhatCC) );
[~,zvalSort] = sort(RhatZ);
rvals = corrcoef([Ytest YnullTest]);
RhatCC_null = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZ_null = 0.5*log( (1+RhatCC_null)./(1-RhatCC_null) );

A = Beta*Vhat;
[~,sortInd] = sort(spref(~isnan(interpolatedR(1,1,1,:))));
sp = spref(~isnan(interpolatedR(1,1,1,:)));
sw = swidth(~isnan(interpolatedR(1,1,1,:)));
spref2 = spref(~isnan(interpolatedR(1,1,1,:)));
spref2 = spref2(sortInd);
A = A(sortInd,:);


RhatScaled = (Rhat.*std(RtoFit,[],[1,2,3]) + mean(RtoFit,[1,2,3]) - mean(R,[1,2,3]))./std(R,[],[1,2,3]);

%% A theoretical treatment, rather than data driven

switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(sp)'-mean(log2(speedsFEF))) ones(size(sp'))];
        
    case 'optimal'
        % More complicated: 'optimal' decoder assuming zero correlations:
        % df/ds*I/(df/ds'*I*df/ds)
        s0 = speedsFEF(speedsFEF==10);
        sprefTemp = sp;
        %     sprefTemp(sprefTemp < 1) = 1;
        swidthTemp = sw;
        %     swidthTemp(swidthTemp > 10) = 10;
        df = -log2(s0./sprefTemp)./(s0.*sw*log(2)).*exp(log2(s0./sprefTemp).^2./(2*sw.^2));
        uOpt = df*inv(eye(length(sprefTemp)))/(df*inv(eye(length(sprefTemp)))*df');
        uOpt(uOpt<-0.05) = min(uOpt(uOpt>-0.05));
        uOptNorm = uOpt/norm(uOpt);
        
        Atheory = [uOpt', ones(size(uOpt'))];
end

% Set up null hypothesis
A0 = randn(size(sp',1),2,1000);

% Normalize 
Atheory = Atheory./vecnorm(Atheory);
A0 = A0./vecnorm(A0);

% Orthogonalize null vector space to theoretical vector space
for nulli = 1:size(A0,3)
    A0(:,1,nulli) = A0(:,1,nulli) - Atheory*inv(Atheory'*Atheory)*Atheory'*A0(:,1,nulli);
    A0(:,2,nulli) = A0(:,2,nulli) - Atheory*inv(Atheory'*Atheory)*Atheory'*A0(:,2,nulli);
end
A0 = A0./vecnorm(A0);

% Construct MT to FEF weight matrix
switch theoretical.expansionDef
    case 'regression'
        XtheoryFit = Xfit*Atheory;
        vOpt = inv(XtheoryFit'*XtheoryFit + eye(size(XtheoryFit'*XtheoryFit)))*XtheoryFit'*Yfit;
        BetaTheory = Atheory*vOpt;
        
        for nulli = 1:size(A0,3)
            X0Fit = Xfit*A0(:,:,nulli);
            v0 = inv(X0Fit'*X0Fit + eye(size(X0Fit'*X0Fit)))*X0Fit'*Yfit;
            Beta0(:,:,nulli) = A0(:,:,nulli)*v0;
        end
        
    case 'bestfit'
        vOpt = inv(Atheory'*Atheory)*Atheory'*BetaRRR;
%         for ai = 1:size(Atheory,2)
%             Wopt(:,:,ai) = Atheory(:,ai)*vOpt(ai,:);
%         end
%         BetaTheory = sum(Wopt,3);
        BetaTheory = Atheory*vOpt;
    case 'data'
        vOpt = Vhat(:,1:size(Atheory,2))';
%         for ai = 1:size(Atheory,2)
%             for ri = 1 %:rankN
%                 tempB(:,:,ri) = Atheory(:,ai)*Vhat(:,ri)';
%             end
%             B(:,:,ai) = sum(tempB,3);
%         end
%         BetaTheory = sum(B,3);
        BetaTheory = Atheory*vOpt;
end
% YtheoryTrain = Xfit*BetaTheory;
YtheoryTest = Xtest*BetaTheory;
Ytheory = Xall*BetaTheory;

for nulli = 1:size(A0,3)
    Y0Test(:,:,nulli) = Xtest*Beta0(:,:,nulli);
end

Rtheory = nan(size(R_z));
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        temp = ismember(conditionsAll,[speedsFEF(si),cohsFEF(ci)],'rows');
        Rtheory(:,si,ci,:) = Ytheory(temp,:);
%         if trainCondition(si,ci)
%             temp = ismember(conditions,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
%             Rtheory(:,si,ci,:) = YtheoryTrain(temp,:);
%         else
%             temp = ismember(conditionsTest,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
%             Rtheory(:,si,ci,:) = YtheoryTest(temp,:);
%         end
    end
end
rvals = corrcoef([Ytest YtheoryTest]);
RhatCCTheory = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZTheory = 0.5*log( (1+RhatCCTheory)./(1-RhatCCTheory) );

for nulli = 1:size(A0,3)
    rvals = corrcoef([Ytest Y0Test(:,:,nulli)]);
    RhatCC0 = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
    RhatZ0(:,nulli) = 0.5*log( (1+RhatCC0)./(1-RhatCC0) );
end

% Compute output along theoretical channels
for ri = 1:size(Atheory,2)
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            theoreticalOutputs(:,si,ci,ri) = inputs_z(:,:,si,ci)'*Atheory(:,ri);
        end
    end
end

RtheoryScaled = (Rtheory.*std(RtoFit,[],[1,2,3]) + mean(RtoFit,[1,2,3]) - mean(R,[1,2,3]))./std(R,[],[1,2,3]);

weightedAverageTheory(:,:,:,1) = sum(R_z./permute(vOpt(1,:),[1,4,3,2]),4);
weightedAverageTheory(:,:,:,2) = sum(R_z./permute(vOpt(2,:),[1,4,3,2]),4);

%% Compare weigth matrix to that expected by theory
q = corrcoef([Atheory(:,1) Beta]);
q = q(1,2:end);
[~,qsort] = sort(q);


%% Analysis of variance
% Signal dependent noise model: Gaussian of the form exp(-w^2/(2*m*ln(s^2))/sqrt(2*pi*m*ln(s^2))
% rankN = 5;
mhat = sqrt(sum(A(spref2 >= 1,1:rankOpt).^2./log(spref2(spref2>=1)'.^2),[1,2])./length(spref2(spref2>=1))); % ML estimator of m 
LLsignal_dependent_noise = sum( -log( sqrt(2*pi*log(spref2(spref2>=1)'.^2)) ) -...
    log(mhat) - A(spref2>=1,1:rankOpt).^2./(2*mhat^2*log(spref2(spref2>=1)'.^2)) ,[1,2]);


% Standard Gaussian noise
sigWeights = sqrt(mean(A(spref2>=1,1:rankOpt).^2,[1,2]));
LLstandard_noise = sum( -log( sqrt(2*pi) ) -...
    log(sigWeights) - A(spref2>=1,1:rankOpt).^2./(2*sigWeights^2) ,[1,2]);

%% Calculate gain from behavior
if compToBehavioralGain.On
    temp = load(compToBehavioralGain.file,'init');
    eye_t = temp.init.t;
    eyeSpeed = temp.init.eye.mean;
    slips = -squeeze(temp.init.eye.mean(temp.init.t == 750,:,:)) + speedsFEF';
    initGain = (temp.init.eye.pert.res - temp.init.eye.pert.resControl)./(0.4*repmat(speedsFEF',[1,3,3]));
    if compToBehavioralGain.applyLogrithmicEstimatorCorrection
        initGain(:,:,2) = initGain(:,:,2).*speedsFEF'*0.4./log2(1.4);
        initGain(:,:,3) = initGain(:,:,3).*speedsFEF'*0.4./log2(1+0.4.*speedsFEF'./slips);
    end
end

%% Saving
if saveResults
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/' dcp{1}.sname ...
        '/MTtoFEFregressionResults'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/MTtoFEFregressionResults' datestr(now,'yyyymmdd')],'-v7.3')
end

%% Plotting
if plotOpts.On

    
    %% FEFsem firing rates
    hFEFrates = figure('Name','FEFsem firing rates during analysis window');
    ylims = [Inf,-Inf];
    for si = 1:length(speedsFEF)
        for ci = 1:length(cohsFEF)
            subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
            plot(fef_t(fef_t >= tWin(1) & fef_t <= tWin(2)),...
                1000*squeeze(Rinit(fef_t >= tWin(1) & fef_t <= tWin(2),si,ci,:)-mean(Rinit(fef_t<=0,:,:,:),[1,2,3])),...
                'Color',speedColors(si,:))
            ylimTemp = ylim;
            ylims(1) = min([ylims(1),ylimTemp(1)]);
            ylims(2) = max([ylims(2),ylimTemp(2)]);
        end
    end
    for si = 1:length(speedsFEF)
        for ci = 1:length(cohsFEF)
            subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
            ylim(ylims);
            xlabel('Time from motion onset (ms)')
            ylabel('\delta Spikes/s (relative to baseline rate)')
        end
    end
    
%     figure('Name','FEFsem firing rates, heatmap')
%     for si = 1:length(speedsFEF)
%         for ci = 1:length(cohsFEF)
%             subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
%             imagesc(fef_t(fef_t >= tWin(1) & fef_t <= tWin(2)),1:size(Rinit,4),squeeze(Rinit(fef_t >= tWin(1) & fef_t <= tWin(2),si,ci,:))')
%         end
%     end
    
    %% Reduced rank regression analysis
    hPredictionVsActual = figure('Name','Predicted vs actual FEFsem responses');
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(cohsFEF),length(cohsFEF),si+(ci-1)*length(cohsFEF))
            plot(squeeze(R_z(taccept,si,ci,:)),squeeze(Rhat(taccept,si,ci,:)),'o','Color',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si)) ' deg/s, Coh = ' num2str(cohsFEF(ci)) '%'])
            hold on
            plotUnity;
            xlabel('FEF data (z-score)')
            ylabel('RRR prediction (z-score)')
            set(gca,'TickDir','out')
        end
    end
    
    %% Z-transformed r-values
    hPredictionZvals = figure('Name','Correlation coefficient between prediction and FEFsem data');
    plot(RhatZ(zvalSort),'bo','DisplayName','Reduced rank model')
    hold on
    plot(RhatZ_full(zvalSort),'ko','DisplayName','Full rank model')
    plot(RhatZTheory(zvalSort),'ro','DisplayName','Theory-based model')
    plot((1:length(zvalSort))',mean(RhatZ0(zvalSort,:),2)+std(RhatZ0(zvalSort,:),[],2)*1.96,'--','Color',[0.6 0.6 0.6])
    plot((1:length(zvalSort))',mean(RhatZ0(zvalSort,:),2)-std(RhatZ0(zvalSort,:),[],2)*1.96,'--','Color',[0.6 0.6 0.6])
%     plotHorizontal(2.34/sqrt(size(Ytest,1)-size(Xfit,2)));       % Approximate pval of 0.01
    xlabel('FEF neuron #')
    ylabel('Fisher transformed r-values')
    set(gca,'TickDir','out')
    
    hNullDistribution = figure;
    subplot(2,1,1)
    histH = histogram(RhatZ'*RhatZ0,200);
    hold on
    theoryH = plotVertical(RhatZ'*RhatZTheory);
    trueH = plotVertical(RhatZ'*RhatZ);
    legend([histH,theoryH,trueH],{'Null','Theory','Actual'})
    xlabel('z * z_{theory}')
    ylabel('Frequency')
    
    subplot(2,1,2)
    histH = histogram(RhatZ'*RhatZ0./vecnorm(RhatZ0)/vecnorm(RhatZ),200);
    hold on
    theoryH = plotVertical(RhatZ'*RhatZTheory/vecnorm(RhatZTheory)/vecnorm(RhatZ));
    legend([histH,theoryH],{'Null','Theory'})
    xlabel('Overlap with reduced rank fit results')
    ylabel('Frequency')
    
    %% Estimate of rank of MT to FEF functional connnections
    hOverallRankCC = figure('Name','Overall rank of MT to FEF weight matrix');
    plot(fitRank(1:10,1),'ko-','DisplayName','Training data')
    hold on
    plot(fitRank(1:10,2),'bo-','DisplayName','Test data')
    plotVertical(rankOpt);
    xlabel('Rank of MT to FEF weight matrix')
    ylabel('Performance')
    set(gca,'TickDir','out')
        
    %%
    hMTinputChannels = figure('Name','MT input channels found by RRR','Position',[674 218 1837 1104]);
    if rankOpt < 5
        rankMax = rankOpt;
    else
        rankMax = 5;
    end
    for ri = 1:rankMax
        for ci = 1:length(cohsFEF)
            subplot(rankMax,length(cohsFEF),ci + (ri-1)*length(cohsFEF))
            for si = 1:length(speedsFEF)
                MToutputChan(:,si,ci,ri) = inputs_z(:,:,si,ci)'*Beta*Vhat(:,ri);
                plot(t,MToutputChan(:,si,ci,ri),'Color',speedColors(si,:),...
                    'DisplayName',['Speed = ' num2str(speedsFEF(si)) ' deg/s, Coh = ' num2str(cohsFEF(ci)) '%'])
                hold on
            end
            xlabel('Time from motion onset (ms)')
            ylabel(['Input along dimension#' num2str(ri)])
            title(['Coh = ' num2str(cohsFEF(ci)) '%'])
            tempLims(ci,:) = ylim;
        end
        
        for ci = 1:length(cohsFEF)
            subplot(rankMax,length(cohsFEF),ci + (ri-1)*length(cohsFEF))
            ylim([min(tempLims(:,1)) max(tempLims(:,2))])
            set(gca,'TickDir','out')
        end
    end
    
    %% Example neurons
    sampleN = 3;
    samplePlots = 1:sampleN*4;
    nonsamplePlots = samplePlots(logical(mod(samplePlots,2)));
    samplePlots = samplePlots(~logical(mod(samplePlots,2)));
    [~,neuronInds] = sort(zvalSort);
    
    hExampleFits = figure('Name','Example fits','Position',[594 119 1917 1203]);
    subplot(sampleN*2,2,nonsamplePlots)
    plot(neuronInds,RhatZ,'bo','DisplayName','Reduced rank model')
    hold on
    xlabel('FEF neuron #')
    ylabel('Fisher transformed r-values')
    set(gca,'TickDir','out')
        
    sampleIndsTemp = randsample(size(Rhat,4)-30:size(Rhat,4),sampleN);
    for sampi = 1:sampleN
        sampleInds(sampi) = find(neuronInds == sampleIndsTemp(sampi));
        plot(neuronInds(sampleInds(sampi)),RhatZ(sampleInds(sampi)),...
            'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0])
        text(neuronInds(sampleInds(sampi)),RhatZ(sampleInds(sampi)),...
            ['Neuron ' num2str(cellID(sampleInds(sampi),1,1:2))])
    end
    
    for sampi = 1:sampleN
        subplot(sampleN*2,2,samplePlots(sampi))
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                if ~trainCondition(si,ci)
                    plot(t(taccept),RtoFit_z(taccept,si,ci,sampleInds(sampi)),...
                        'Color',[speedColors(si,:) 1-cohsFEF(ci)/100],...
                        'DisplayName',['Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
                    hold on
                    plot(t(taccept),Rhat(taccept,si,ci,sampleInds(sampi)),...
                        'k--',...
                        'DisplayName',['Model, Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
                end
            end
        end
        title(['Neuron ' num2str(cellID(sampleInds(sampi),1,1:2))])
        xlabel('Time from motion onset (ms)')
        ylabel('z-score')
        set(gca,'TickDir','out')
    end
    
    
    subplot(sampleN*2,2,nonsamplePlots)
    sampleIndsTemp = randsample(1:30,sampleN);
    for sampi = 1:sampleN
        sampleInds(sampi) = find(neuronInds == sampleIndsTemp(sampi));
        plot(neuronInds(sampleInds(sampi)),RhatZ(sampleInds(sampi)),...
            'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0])
        text(neuronInds(sampleInds(sampi)),RhatZ(sampleInds(sampi)),...
            ['Neuron ' num2str(cellID(sampleInds(sampi),1,1:2))])
    end
    
    for sampi = 1:sampleN
        subplot(sampleN*2,2,samplePlots(sampi)+2*sampleN)
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                if ~trainCondition(si,ci)
                    plot(t(taccept),RtoFit_z(taccept,si,ci,sampleInds(sampi)),...
                        'Color',[speedColors(si,:) 1-cohsFEF(ci)/100],...
                        'DisplayName',['Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
                    hold on
                    plot(t(taccept),Rhat(taccept,si,ci,sampleInds(sampi)),...
                        'k--',...
                        'DisplayName',['Model, Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
                end
            end
        end
        title(['Neuron ' num2str(cellID(sampleInds(sampi),1,1:2))])
        xlabel('Time from motion onset (ms)')
        ylabel('z-score')
        set(gca,'TickDir','out')
    end
    
    %% Mean and mean fit values
    hMeanFitInitiation = figure('Name','Comparison of means','Position',[391 274 1778 990]);
    lims = [Inf,-Inf];
    for ci = 1:length(cohsFEF)
        subplot(2,length(cohsFEF),ci)
        for si = 1:length(speedsFEF)
            plot(t(taccept),mean(RtoFit_z(:,si,ci,:),4),'Color',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
            hold on
            plot(t(taccept),mean(Rhat(taccept,si,ci,:),4),'k--',...
                'DisplayName',['Reduced rank model; Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
        end
        axis tight
        templim = ylim;
        lims(1) = min([lims(1),templim(1)]);
        lims(2) = max([lims(2),templim(2)]);
        xlabel('Time from motion onset (ms)')
        ylabel('Mean z-scored response')
        title(['Coh = ' num2str(cohsFEF(ci)) '%'])
    end
    for ci = 1:length(cohsFEF)
        subplot(2,length(cohsFEF),length(cohsFEF)+ci)
        for si = 1:length(speedsFEF)
            plot(t(taccept),mean(RtoFit_z(:,si,ci,:),4),'Color',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
            hold on
            plot(t(taccept),mean(Rtheory(taccept,si,ci,:),4),'r--',...
                'DisplayName',['Theory-based model; Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
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
        set(gca,'TickDir','out')
    end
    
    hMeanFit = figure('Name','Comparison of means','Position',[391 274 1778 990]);
    for ci = 1:length(cohsFEF)
        subplot(2,length(cohsFEF),ci)
        for si = 1:length(speedsFEF)
            plot(t,mean(R_z(:,si,ci,:),4),'Color',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
            hold on
            plot(t,mean(RhatScaled(:,si,ci,:),4),'k--',...
                'DisplayName',['Reduced rank model; Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
        end
        axis tight
        templim = ylim;
        lims(1) = min([lims(1),templim(1)]);
        lims(2) = max([lims(2),templim(2)]);
        xlabel('Time from motion onset (ms)')
        ylabel('Mean z-scored response')
        title(['Coh = ' num2str(cohsFEF(ci)) '%'])
    end
    
    for ci = 1:length(cohsFEF)
        subplot(2,length(cohsFEF),length(cohsFEF)+ci)
        for si = 1:length(speedsFEF)
            plot(t,mean(R_z(:,si,ci,:),4),'Color',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
            hold on
            plot(t,mean(RtheoryScaled(:,si,ci,:),4),'r--',...
                'DisplayName',['Theory-based model; Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci)) '%'])
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
        set(gca,'TickDir','out')
    end
    
    %% Weights along 1st reduced rank dimension as a function of preferred speed
    figure
    ranks = 1:rankOpt;
    for ri = 1:length(ranks)
        subplot(length(ranks),2,1+(ri-1)*2)
        semilogx(spref2(spref2>=1),A(spref2>=1,ranks(ri)),'o')
        hold on
        semilogx(sp(sp>=1),Atheory(sp>=1,1),'o')
        s = logspace(log10(speedPrefOpts.lb(1)),log10(speedPrefOpts.ub(1)),30);
        if LLsignal_dependent_noise > LLstandard_noise
            plot(s,mhat*sqrt(log(s.^2)),'r')
            plot(s,-mhat*sqrt(log(s.^2)),'r')
        end
        
        plot(s,sigWeights*ones(length(s),1),'k')
        plot(s,-sigWeights*ones(length(s),1),'k')
        
        xlim([speedPrefOpts.lb(1), speedPrefOpts.ub(1)])
        ax = axis;
        text(2,0.95*ax(4),['\delta LL = ' num2str(LLsignal_dependent_noise-LLstandard_noise)])
        
        xlabel('Preferred speed (deg/s)')
        ylabel('Weight')
        
        subplot(length(ranks),2,2+(ri-1)*2)
        plot(Vhat(zvalSort,ranks(ri)),'o')
    end
    
    %% Compare output channels to behavioral gain
    if compToBehavioralGain.On
        
        figure('Name','MT output vs. gain')
        if rankOpt < 5
            rankMax = rankOpt;
        else
            rankMax = 5;
        end
        for ri = 1:rankMax
                subplot(1,rankMax,ri)
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
        
        
        figure('Name','Theoretical MT output vs. gain')
        for ri = 1:size(Atheory,2)
                subplot(1,size(Atheory,2),ri)
                for si = 1:length(speedsFEF)                    
                    outputTemp = squeeze(nanmean(theoreticalOutputs(t >= 700 & t <= 800,si,:,ri),1));
                    plot(squeeze(initGain(si,:,3)),outputTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                    hold on
                    outputTemp = squeeze(nanmean(theoreticalOutputs(t >= 100 & t <= 200,si,:,ri),1));
                    plot(squeeze(initGain(si,:,2)),outputTemp,...
                        'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                end
                xlabel('Behavioral gain (unitless)')
                ylabel(['MT input along dimension ' num2str(ri)])
        end
        
        hFEFvsTargeSpeedEyeSpeedGain = figure('Name','FEF output vs target speed, eye speed and gain','Position',[543 432 1646 413]);
        subplot(1,3,1)
        for ci = 1:length(cohsFEF)
            if tWin(2) < 200
                outputTemp = squeeze(nanmean(Rinit(t >= tWin(2) & t <= tWin(2),:,ci,:),[1,4]));
            else
                outputTemp = squeeze(nanmean(Rinit(t >= 100 & t <= 200,:,ci,:),[1,4]));
            end
            plot(speedsFEF,outputTemp*1000,...
                'd-','Color',1-cohsFEF(ci)/100*ones(1,3),'MarkerFaceColor',1-cohsFEF(ci)/100*ones(1,3))
            hold on
        end
        xlabel('Target speed (deg/s)')
        ylabel(['Mean FEFsem activity (sp/s)'])
        axis square
        set(gca,'TickDir','out')
        
        subplot(1,3,2)
        for si = 1:length(speedsFEF)
            if tWin(2) < 200
                outputTemp = squeeze(nanmean(Rinit(t >= tWin(2) & t <= tWin(2),si,:,:),[1,4]));
            else
                outputTemp = squeeze(nanmean(Rinit(t >= 100 & t <= 200,si,:,:),[1,4]));
            end
            plot(squeeze(eyeSpeed(eye_t == tWin(2),si,:)),outputTemp*1000,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            hold on
        end
        xlabel('Eye speed (deg/s)')
        ylabel(['Mean FEFsem activity (sp/s)'])
        axis square
        set(gca,'TickDir','out')
        
        subplot(1,3,3)
        for si = 1:length(speedsFEF)
            if tWin(2) < 200
                outputTemp = squeeze(nanmean(Rinit(t >= tWin(2) & t <= tWin(2),si,:,:),[1,4]));
            else
                outputTemp = squeeze(nanmean(Rinit(t >= 100 & t <= 200,si,:,:),[1,4]));
            end
            plot(squeeze(initGain(si,:,2)),outputTemp*1000,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            hold on
        end
        xlabel('Behavioral gain (unitless)')
        ylabel(['Mean FEFsem activity (sp/s)'])
        axis square
        set(gca,'TickDir','out')
    end
    
    %% Plot reconstructed firing rates projection on gain dimension
%     if compToBehavioralGain.On
%         temp = load(compToBehavioralGain.file,'BinitOrth');
%         BinitOrth = temp.BinitOrth;
%         sz = size(Rhat);
%         pInit = reshape((BinitOrth(:,1:4)'*reshape(Rhat,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
%         pInitTheory = reshape((BinitOrth(:,1:4)'*reshape(Rtheory,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
%         
%         figure('Name','Predicted response along targeted dimensions','Position',[63 169 1606 1079])
%         subplot(3,1,1)
%         for speedi = 1:length(speedsFEF)
%             for cohi = 1:length(cohsFEF)
%                 plot(t,pInit(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%         end
%         axis tight
%         ax(1,:) = axis;
%         xlabel('Time from motion onset (ms)')
%         ylabel('Speed related activity (a.u.)')
%         
%         subplot(3,1,2)
%         for cohi = 1:length(cohsFEF)
%             for speedi = 1:length(speedsFEF)
%                 plot(t,pInit(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%         end
%         axis tight
%         ax(1,:) = axis;
%         xlabel('Time from motion onset (ms)')
%         ylabel('Coherence related activity (a.u.)')
%         
%         subplot(3,1,3)
%         for cohi = 1:length(cohsFEF)
%             for speedi = 1:length(speedsFEF)
%                 plot(t,pInit(:,speedi,cohi,3),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%         end
%         axis tight
%         ax(1,:) = axis;
%         xlabel('Time from motion onset (ms)')
%         ylabel('Gain related activity (a.u.)')
%         
%         
%         dimNames = {'Speed','Coherence','Gain','Offset'};
%         gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
%         for targetedDim = 1:3
%             tempData = [];
%             subplot(3,1,targetedDim)
%             for si = 1:length(speedsFEF)
%                 initRatesTemp = squeeze(nanmean(pInit(t >= 700 & t <= 800,si,:,targetedDim),1));
%                 plot(squeeze(initGain(si,:,3)),initRatesTemp,...
%                     'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%                 hold on
%                 tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%                 initRatesTemp = squeeze(nanmean(pInit(t >= 100 & t <= 200,si,:,targetedDim),1));
%                 plot(squeeze(initGain(si,:,2)),initRatesTemp,...
%                     'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%                 tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
%             end
%             gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
%             axis square
%             ax = axis;
%             text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
%             xlabel('Behavioral gain (unitless)')
%             ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
%         end
%         
%         figure('Name','Theory-based response along targeted dimensions','Position',[63 169 1606 1079])
%         subplot(3,1,1)
%         for speedi = 1:length(speedsFEF)
%             for cohi = 1:length(cohsFEF)
%                 plot(t,pInitTheory(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%         end
%         axis tight
%         ax(1,:) = axis;
%         xlabel('Time from motion onset (ms)')
%         ylabel('Speed related activity (a.u.)')
%         
%         subplot(3,1,2)
%         for cohi = 1:length(cohsFEF)
%             for speedi = 1:length(speedsFEF)
%                 plot(t,pInitTheory(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%         end
%         axis tight
%         ax(1,:) = axis;
%         xlabel('Time from motion onset (ms)')
%         ylabel('Coherence related activity (a.u.)')
%         
%         subplot(3,1,3)
%         for cohi = 1:length(cohsFEF)
%             for speedi = 1:length(speedsFEF)
%                 plot(t,pInitTheory(:,speedi,cohi,3),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%         end
%         axis tight
%         ax(1,:) = axis;
%         xlabel('Time from motion onset (ms)')
%         ylabel('Gain related activity (a.u.)')
%         
%         
%         dimNames = {'Speed','Coherence','Gain','Offset'};
%         gvth2 = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
%         for targetedDim = 1:3
%             tempData = [];
%             subplot(3,1,targetedDim)
%             for si = 1:length(speedsFEF)
%                 initRatesTemp = squeeze(nanmean(pInitTheory(t >= 700 & t <= 800,si,:,targetedDim),1));
%                 plot(squeeze(initGain(si,:,3)),initRatesTemp,...
%                     'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%                 hold on
%                 tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%                 initRatesTemp = squeeze(nanmean(pInitTheory(t >= 100 & t <= 200,si,:,targetedDim),1));
%                 plot(squeeze(initGain(si,:,2)),initRatesTemp,...
%                     'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%                 tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
%             end
%             gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
%             axis square
%             ax = axis;
%             text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
%             xlabel('Behavioral gain (unitless)')
%             ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
%         end
%     end
    
    %% Relationship to Functional topography
    if neuronTypingComparison.On
        temp = load(neuronTypingComparison.file,'NumClusters','Y','idx','centInd','cellID','subjects');
        subjectInd = find(strcmp(temp.subjects,dcp{1}.sname));
        Y = temp.Y;
        Y = Y(temp.cellID(:,1,4) == subjectInd,:);
        idx = temp.idx(temp.cellID(:,1,4)==subjectInd);
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
                
      
%         figure
%         ranks = 1:4;
%         neuronN = 1:size(Vhat,1);
%         for ri = 1:length(ranks)
%             subplot(1,length(ranks),ri)
%             for typei = 1:NumClusters
%                 plot(neuronN(temp.idx == typei),Vhat(temp.idx == typei,ri),'o',...
%                     'Color',colorWheel(typei,:),'MarkerFaceColor',colorWheel(typei,:))
%                 hold on
%             end
%             xlabel('FEF neuron number')
%             ylabel(['Loading in dim ' num2str(ri)])
%         end
        
        hFitQualityVsFunctionalTopology = figure('Name','Fit quality vs functional topology','Position',[599 758 1392 413]);
        subplot(1,2,1)
        scatter(Y(:,1),Y(:,2),abs(RhatZ)*100,RhatZ,'filled')
        cax(1,:) = caxis;
        axis square
        xlabel('tSNE_1')
        ylabel('tSNE_2')
        
        subplot(1,2,2)
        scatter(Y(:,1),Y(:,2),abs(RhatZTheory)*100,RhatZTheory,'filled')
        axis square
        cax(2,:) = caxis;
        xlabel('tSNE_1')
        ylabel('tSNE_2')
        
        subplot(1,2,1)
        caxis([min(cax(:,1)) max(cax(:,2))])
        set(gca,'TickDir','out')
        subplot(1,2,2)
        caxis([min(cax(:,1)) max(cax(:,2))])
        set(gca,'TickDir','out')
        
        figure
        idxTemp = idx(qsort);
        for typei = 1:NumClusters
            subplot(3,NumClusters,typei)
            for si = 1:length(speedsFEF)
                for ci = 1:length(cohsFEF)
                    plot(t,nanmean(R_z(:,si,ci,idx==typei),4)*1000,...
                        'Color',[speedColors(si,:) cohsFEF(ci)/100])
                    hold on
                end
            end
            xlabel('Time from motion onset (ms)')
            ylabel('sp/s')
            
            subplot(3,NumClusters,NumClusters+typei)
            imagesc(BetaRRR(sortInd,idxTemp==typei))
            ylabel('MT unit (sorted by speed pref)')
            xlabel('FEF unit in cluster')
            
            subplot(3,NumClusters,2*NumClusters+typei)
            Atemp = Beta(:,idxTemp==typei)*Vhat(idxTemp==typei,:);
            semilogx(spref2,Atemp(sortInd,1),'o')
            ylabel('MT unit weight')
            xlabel('speed preference (deg/s)')
        end
        
        
        figure
        for typei = 1:NumClusters
            subplot(3,NumClusters,typei)
            for si = 1:length(speedsFEF)
                for ci = 1:length(cohsFEF)
                    plot(t,nanmean(Rtheory(:,si,ci,idx==typei),4)*1000,...
                        'Color',[speedColors(si,:) cohsFEF(ci)/100])
                    hold on
                end
            end
            xlabel('Time from motion onset (ms)')
            ylabel('sp/s')
            
            subplot(3,NumClusters,NumClusters+typei)
            imagesc(BetaTheory(sortInd,idxTemp==typei))
            ylabel('MT unit (sorted by speed pref)')
            xlabel('FEF unit in cluster')
            
            subplot(3,NumClusters,2*NumClusters+typei)
            Atemp = BetaTheory(:,idxTemp==typei)*Vhat(idxTemp==typei,:);
            semilogx(spref2,Atemp(sortInd,1:size(Atheory,2)),'o')
            ylabel('MT unit weight')
            xlabel('speed preference (deg/s)')
        end
    end
    
    %% Theory guided average
%     if compToBehavioralGain.On
%         
%         figure('Name','Predicted response, weighted average','Position',[73 206 1606 700])
%         for speedi = 1:length(speedsFEF)
%             subplot(3,3,speedi)
%             for cohi = 1:length(cohsFEF)
%                 plot(t,weightedAverageTheory(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%             axis tight
%             ax(speedi,:) = axis;
%             xlabel('Time from motion onset (ms)')
%             ylabel('Mean activity, weighred by theoretical results (a.u.)')
%         end
%         for speedi = 1:length(speedsFEF)
%             subplot(3,3,speedi)
%             ylim([min(ax(:,3)) max(ax(:,4))])
%             plotVertical(tempWin);
%         end
%         for speedi = 1:length(speedsFEF)
%             subplot(3,3,speedi+length(speedsFEF))
%             for cohi = 1:length(cohsFEF)
%                 plot(t,weightedAverageTheory(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%             axis tight
%             ax(speedi,:) = axis;
%             xlabel('Time from motion onset (ms)')
%             ylabel('Mean activity, weighred by theoretical results (a.u.)')
%         end
%         for speedi = 1:length(speedsFEF)
%             subplot(3,3,speedi+length(speedsFEF))
%             ylim([min(ax(:,3)) max(ax(:,4))])
%             plotVertical(tempWin);
%         end
%         for speedi = 1:length(speedsFEF)
%             subplot(3,3,speedi+2*length(speedsFEF))
%             for cohi = 1:length(cohsFEF)
%                 plot(t,sum(weightedAverageTheory(:,speedi,cohi,:),4),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
%                 hold on
%             end
%             axis tight
%             ax(speedi,:) = axis;
%             xlabel('Time from motion onset (ms)')
%             ylabel('Mean activity, weighred by theoretical results (a.u.)')
%         end
%         for speedi = 1:length(speedsFEF)
%             subplot(3,3,speedi+2*length(speedsFEF))
%             ylim([min(ax(:,3)) max(ax(:,4))])
%             plotVertical(tempWin);
%         end
%         
%         gvth3 = figure('Name',['Behavioral gain vs activity weighted by their relationships to the theoretical gain computation'],'Position',[1956 59 570 1263]);
%         tempData = [];        
%         subplot(2,1,1)
%         for si = 1:length(speedsFEF)
%             initRatesTemp = squeeze(nanmean(weightedAverageTheory(t >= 700 & t <= 800,si,:,:),[1,4]));
%             plot(squeeze(initGain(si,:,3)),initRatesTemp,...
%                 'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%             hold on
%             tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%             initRatesTemp = squeeze(nanmean(weightedAverageTheory(t >= tempWin(1) & t <= tempWin(2),si,:,:),[1,4]));
%             plot(squeeze(initGain(si,:,2)),initRatesTemp,...
%                 'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%             tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
%         end
%         gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
%         axis square
%         ax = axis;
%         text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
%         xlabel('Behavioral gain (unitless)')
%         ylabel(['Mean of activity, weighted by theoretical results'])
%         
%         tempData = [];
%         subplot(2,1,2)
%         for si = 1:length(speedsFEF)
%             initRatesTemp = squeeze(nanmean(weightedAverageTheory(t >= 700 & t <= 800,si,:,:),[1,4]));
%             plot(squeeze(initGain(si,:,3)),initRatesTemp,...
%                 'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%             hold on
%             tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
%             initRatesTemp = squeeze(nanmean(weightedAverageTheory(t >= tempWin(1) & t <= tempWin(2),si,:,:),[1,4]));
%             plot(squeeze(initGain(si,:,2)),initRatesTemp,...
%                 'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
%             tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
%         end
%         gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
%         axis square
%         ax = axis;
%         text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
%         xlabel('Behavioral gain (corrected, unitless)')
%         ylabel(['Mean of activity, weighted by theoretical results'])
%         
%     end
    %% Save figures
    if saveFigures

        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' dcp{1}.sname ...
            '/MTtoFEFregression/' datestr(now,'yyyymmdd')];
        if ~exist(saveLocation,'dir')
            mkdir(saveLocation)
        end

        savefig(hFEFrates,[saveLocation '/FEFfiringInWindow.fig'])
        savefig(hPredictionVsActual,[saveLocation '/PredictionVsActual.fig'])
        savefig(hPredictionZvals ,[saveLocation '/ZtransformedCorrelationsByNeuron.fig'])
        savefig(hOverallRankCC ,[saveLocation '/RankOfMTtoFEFmatrix.fig'])
        savefig(hMTinputChannels ,[saveLocation '/InputsAlongMTchannels.fig'])
        savefig(hMeanFitInitiation ,[saveLocation '/meanFEFinitiationVsModels.fig'])
        savefig(hExampleFits,[saveLocation '/exampleFits.fig'])

        if compToBehavioralGain.On
            savefig(hFEFvsTargeSpeedEyeSpeedGain,[saveLocation '/FEFvsTargetSpeedEyeSpeedGain.fig'])
        end

        if neuronTypingComparison.On
            savefig(hFitQualityVsFunctionalTopology,[saveLocation '/MTtoFEFfitQualityFunctionalTopology.fig'])
        end

    end

end

