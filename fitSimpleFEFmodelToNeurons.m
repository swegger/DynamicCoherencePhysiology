function modelFEF = fitSimpleFEFmodelToNeurons(dcp,varargin)
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
addParameter(Parser,'tWin',[-100 900])
addParameter(Parser,'tau',20)
addParameter(Parser,'lambdaRidge',0)
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'weightPrior',NaN)
addParameter(Parser,'tauLB',NaN)
addParameter(Parser,'tauUB',NaN)
addParameter(Parser,'theoretical',theoretical_default)
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
tWin = Parser.Results.tWin;
tau = Parser.Results.tau;
lambdaRidge = Parser.Results.lambdaRidge;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
weightPrior = Parser.Results.weightPrior;
tauLB = Parser.Results.tauLB;
tauUB = Parser.Results.tauUB;
theoretical = Parser.Results.theoretical;
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

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff);
cellID = cellID(passCutoff,:,:);

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
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);
spref2 = spref(prod(any(isnan(inputs),2),[3,4])==0);
swidth2 = swidth(prod(any(isnan(inputs),2),[3,4])==0);

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
    P0 = [tauLB/dt,0,(log2(spref2)-mean(log2(speedsFEF)))/1000];
end
if any(isnan(lb))
    lb = [tauLB/dt, -100, -Inf*ones(1,size(inputs,1))];
end
if any(isnan(ub))
    ub = [tauUB/dt, 100, Inf*ones(1,size(inputs,1))];
end

taccept = t >= tWin(1) & t <= tWin(2);
mask = ones(size(trainCondition));
mask(~trainCondition) = NaN;

% OPTIONS = optimset(@fmincon);
% OPTIONS.MaxFunEvals = 3e10;
% OPTIONS.MaxIterations = 1e10;
OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10);
for neuroni = 1:size(R,4)
    tic
    disp(['Neuron ' num2str(neuroni) ' of ' num2str(size(R,4))])
    Rtemp = R_z(taccept,:,:,neuroni);
    Rtemp = Rtemp.*permute(mask,[3,1,2]);
    inputsTemp = inputs_z(:,taccept,:,:);
    inputsTemp = inputsTemp.*permute(mask,[3,4,1,2]);
    R0 = nanmean(Rtemp(1,:,:),[2,3]);
    minimizer = @(P)minimizant(P,R0,inputsTemp,Rtemp,lambdaRidge,W0);
    [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
    
    modelFEF.tau(neuroni) = P(1)*dt;
    modelFEF.R0(neuroni) = R0;
    modelFEF.baseLine(neuroni) = P(2);
    modelFEF.W(:,neuroni) = P(3:end);
    modelFEF.fval(:,neuroni) = fval;
    modelFEF.exitflag(:,neuroni) = exitflag;
    modelFEF.output(:,neuroni) = output;
    modelFEF.lambda(:,neuroni) = lambda;
    modelFEF.grad(:,neuroni) = grad;
    modelFEF.hessian(:,:,neuroni) = hessian;
    toc
end

modelFEF.dt = dt;


%% Analysis

% Sort MT dimension by speed preference, remove poor fit FEF units, and
% remove leak rate
[~,sortInd] = sort(spref(~isnan(interpolatedR(1,1,1,:))));
spref2 = spref(~isnan(interpolatedR(1,1,1,:)));
spref2 = spref2(sortInd);
cutThreshold = quantile(modelFEF.fval,0.975);
W = modelFEF.W;

% Normalize by maximium value of W for each FEF neuron
wmax = max(W,[],2);
Weff = W./repmat(wmax,[1,size(W,2)]);

% SVD
[U,S,V] = svd(Weff);

% Build matrices
X = nan(size(inputs,2)*sum(trainCondition(:)),size(inputs,1));
Y = nan(size(R_z,1)*sum(trainCondition(:)),size(R_z,4));
conditions = nan(size(R_z,1)*sum(trainCondition(:)),2);

idx = [1, 1];
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        X(idx(1):idx(1)+size(inputs,2)-1,:) = permute(inputs_z(:,:,si,ci),[2,1]);
        Y(idx(2):idx(2)+size(R_z,1)-1,:) = squeeze(R_z(:,si,ci,:));
        conditions(idx(2):idx(2)+size(R_z,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci)],[size(R_z,1),1]);
        idx(1) = idx(1) + size(inputs,2);
        idx(2) = idx(2) + size(R_z,1);
    end
end
[Uhat,Shat,Vhat] = svd(Y,0);
A = W*Vhat;

%% A theoretical treatment, rather than data driven

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
        df = -log2(s0./sprefTemp)./(s0.*sw*log(2)).*exp(log2(s0./sprefTemp).^2./(2*sw.^2));
        uOpt = df*inv(eye(length(sprefTemp)))/(df*inv(eye(length(sprefTemp)))*df');
        uOpt(uOpt<-0.05) = min(uOpt(uOpt>-0.05));
        uOptNorm = uOpt/norm(uOpt);
        
        Atheory = [uOpt', ones(size(uOpt'))];
end

% Normalize 
Atheory = Atheory./vecnorm(Atheory);

% Construct MT to FEF weight matrix
switch theoretical.expansionDef
    case 'bestfit'
        vOpt = inv(Atheory'*Atheory)*Atheory'*W;
        for ai = 1:size(Atheory,2)
            Wopt(:,:,ai) = Atheory(:,ai)*vOpt(ai,:);
        end
        WTheory = sum(Wopt,3);
    case 'data'
        for ai = 1:size(Atheory,2)
            for ri = 1 %:rankN
                tempB(:,:,ri) = Atheory(:,ai)*Vhat(:,ri)';
            end
            B(:,:,ai) = sum(tempB,3);
        end
        WTheory = sum(B,3);
end

%% Make predictions
for neuroni = size(R_z,4)
    Rhat(:,:,:,neuroni) = SimpleFEFmodel(modelFEF.W(:,neuroni)',modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),modelFEF.tau(neuroni)/modelFEF.dt,inputs);
    RhatTheory(:,:,:,neuroni) = SimpleFEFmodel(WTheory(:,neuroni)',0,modelFEF.R0(neuroni),modelFEF.tau(neuroni)/modelFEF.dt,inputs_z);
end


%% Plotting
if plotOpts.On

    %% Fit quality
    figure
    for neuroni = size(R,4)
        tempFEF = R(:,:,:,neuroni);
        test = SimpleFEFmodel(modelFEF.W(:,neuroni)',modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),modelFEF.tau(neuroni)/modelFEF.dt,inputs);
        for ci = 1:length(cohsFEF)
            subplot(1,length(cohsFEF),ci)
            for si = 1:length(speedsFEF)
                plot(tempFEF(taccept,si,ci)*1000,test(taccept,si,ci)*1000,'o','Color',speedColors(si,:))
                hold on
            end
        end
    end
    
    for ci = 1:length(cohsFEF)
        subplot(1,length(cohsFEF),ci)
        %     axis(1000*[min(RtoFit(:))-0.1*min(RtoFit(:)), max(RtoFit(:))+0.01*max(RtoFit(:)), ...
        %         min(RtoFit(:))-0.1*min(RtoFit(:)), max(RtoFit(:))+0.01*max(RtoFit(:))])
        axis square
        plotUnity;
        xlabel('FEF data (sp/s)')
        ylabel('Fit data (sp/s)')
    end

end


%% Saving
if saveOpts.On
    save(saveOpts.location);
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
function out = minimizant(P,R0,inputs,R,lambda,W0)
    tau = P(1);
    baseLine = P(2);
    W = P(3:end);
    Rest = SimpleFEFmodel(W,baseLine,R0,tau,inputs);
    out = nansum((R(:) - Rest(:)).^2) + lambda*nansum((W-W0).^2);