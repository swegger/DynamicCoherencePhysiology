function modelFEF = fitReducedRankModel(subject,R,fef_t,eye_t,initGain,meanEyeSpeed,dimNames,varargin)
%% neuronPartialCorrelation
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
theoretical_default.weightTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';
equalizeInputsPriorToStimulusOnset_default.On = false;
simulateMT_default.On = false;

centerData_default.On = false;

loadFromDataFile_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addRequired(Parser,'R')
addRequired(Parser,'fef_t')
addRequired(Parser,'eye_t')
addRequired(Parser,'initGain')
addRequired(Parser,'meanEyeSpeed')
addRequired(Parser,'dimNames')
addParameter(Parser,'fitType','fitSigmas')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'trainCondition',trainCondition_default)
addParameter(Parser,'P0',NaN)
addParameter(Parser,'ub',NaN)
addParameter(Parser,'lb',NaN)
addParameter(Parser,'sigmaLowerBound',0.05)
addParameter(Parser,'tWin',[-100 900])
addParameter(Parser,'tau',20)
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'simulateMT',simulateMT_default)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'weightPrior',NaN)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset_default)
addParameter(Parser,'centerData',centerData_default)
addParameter(Parser,'probeObjectiveFunction',false)
addParameter(Parser,'plotOpts',plotOpts_default)
addParameter(Parser,'saveResults',false)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'loadFromDataFile',loadFromDataFile_default)
addParameter(Parser,'saveLocation',[])
addParameter(Parser,'reducedSave',false)

parse(Parser,subject,R,fef_t,eye_t,initGain,meanEyeSpeed,dimNames,varargin{:})

subject = Parser.Results.subject;
R = Parser.Results.R;
fef_t = Parser.Results.fef_t;
eye_t = Parser.Results.eye_t;
initGain = Parser.Results.initGain;
meanEyeSpeed = Parser.Results.meanEyeSpeed;
dimNames = Parser.Results.dimNames;
fitType = Parser.Results.fitType;
objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
trainCondition = Parser.Results.trainCondition;
P0 = Parser.Results.P0;
ub = Parser.Results.ub;
lb = Parser.Results.lb;
sigmaLowerBound = Parser.Results.sigmaLowerBound;
tWin = Parser.Results.tWin;
tau = Parser.Results.tau;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
simulateMT = Parser.Results.simulateMT;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
theoretical = Parser.Results.theoretical;
equalizeInputsPriorToStimulusOnset = Parser.Results.equalizeInputsPriorToStimulusOnset;
centerData = Parser.Results.centerData;
probeObjectiveFunction = Parser.Results.probeObjectiveFunction;
plotOpts = Parser.Results.plotOpts;
saveResults = Parser.Results.saveResults;
saveFigures = Parser.Results.saveFigures;
loadFromDataFile = Parser.Results.loadFromDataFile;
saveLocation = Parser.Results.saveLocation;
reducedSave = Parser.Results.reducedSave;

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

if loadFromDataFile.On
    saveResultsOld = saveResults;
    saveFiguresOld = saveFigures;
    plotOptsOld = plotOpts;
    load(loadFromDataFile.File)  
    saveResults = saveResultsOld;
    saveFigures = saveFiguresOld;
    plotOpts = plotOptsOld;
else
%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
if simulateMT.On
    t = -500:900;
    speeds = speedsFEF;
    [MT, spref, cohs] = simulateMTdata(mtResults,t,simulateMT.modelN,...
        'speedsMT',speeds,'removeBaseline',simulateMT.removeBaseline,'gaussianApprox',simulateMT.gaussianApprox,...
        'accept',simulateMT.accept);
    swidth = 1.2*ones(size(spref));
else
    [MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
        'speedsMT',speeds,'cohsMT',cohs,'directionsMT',directionsMT,...
        'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
        'speedPrefOpts',speedPrefOpts);
    
    t = mtResults.mt{1}.neuron_t;
end


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);

inputs = permute(interpolatedR,[4,1,2,3]);
if equalizeInputsPriorToStimulusOnset.On
    switch equalizeInputsPriorToStimulusOnset.method
        case 'conditions'
            inputs(:,t <= equalizeInputsPriorToStimulusOnset.latency,:,:) = ...
                repmat(mean(inputs(:,t <= equalizeInputsPriorToStimulusOnset.latency,:,:),[3,4]),[1,1,size(inputs,3),size(inputs,4)]);
        case 'time&conditions'
            inputs(:,t <= equalizeInputsPriorToStimulusOnset.latency,:,:) = ...
                repmat(mean(inputs(:,t <= equalizeInputsPriorToStimulusOnset.latency,:,:),[2,3,4]),[1,sum(t<=equalizeInputsPriorToStimulusOnset.latency),size(inputs,3),size(inputs,4)]);
        otherwise
            error('Equalization method not recognized.')
    end
end
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);
sp = spref(~isnan(interpolatedR(1,1,1,:)));%spref(prod(any(isnan(inputs),2),[3,4])==0);
sw = swidth(~isnan(interpolatedR(1,1,1,:)));%swidth(prod(any(isnan(inputs),2),[3,4])==0);

%% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,t);
R = R(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);
dt = t(2)-t(1);
t = t(iMT);

inputs_z = (inputs-mean(inputs(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(inputs(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);
inputs_n = inputs./max(abs(inputs),[],[2,3,4]);

R_z = (R-mean(R,[1,2,3]))./std(R,[],[1,2,3]);
R_n = R./max(abs(R),[],[1,2,3]);

taccept = t >= tWin(1) & t <= tWin(2);

RtoFit = R_n(:,:,:,~strcmp(dimNames,'Offset'));
RtoFit = RtoFit - mean(RtoFit(t<=0,:,:,:),[1,2,3]);
if centerData.On
    RtoFit = RtoFit - RtoFit(:,centerData.inds(1),centerData.inds(2),:);
    initGain = initGain - initGain(centerData.inds(1),centerData.inds(2),:);
end
RtoFit = RtoFit./max(abs(RtoFit),[],[1,2,3]);
% RtoFit(:,:,:,2) = -RtoFit(:,:,:,2);

%% Set MT weights
switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(sp)' - log2(mean(speedsFEF))) ones(size(sp'))];
%         Atheory = [(log2(sp)') ones(size(sp'))];
        
    case 'optimal'
        % More complicated: 'optimal' decoder assuming zero correlations:
        % df/ds*I/(df/ds'*I*df/ds)
        s0 = speedsFEF(speedsFEF==10);
        sprefTemp = sp;
        swidthTemp = sw;
        df = -log2(s0./sprefTemp)./(s0.*swidthTemp*log(2)).*exp(log2(s0./sprefTemp).^2./(2*swidthTemp.^2));
        uOpt = df*inv(eye(length(sprefTemp)))/(df*inv(eye(length(sprefTemp)))*df');
        uOpt(uOpt<-0.05) = min(uOpt(uOpt>-0.05));
        uOptNorm = uOpt/norm(uOpt);
        
        Atheory = [uOpt', ones(size(uOpt'))];
end

% Normalize 
Atheory = Atheory./vecnorm(Atheory);

% Compute MT output channels
for ri = 1:size(Atheory,2)
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            theoreticalInput(ri,:,si,ci) = inputs_n(:,:,si,ci)'*Atheory(:,ri);
        end
    end
end
if any(strcmp(dimNames,'Offset'))
    theoreticalInput(size(Atheory,2)+1,:,:,:) = repmat(permute(mean(R_n(:,:,:,strcmp(dimNames,'Offset')),[2,3]),[2,1]),[1,1,length(speedsFEF),length(cohsFEF)]);
end
theoreticalInput = theoreticalInput - theoreticalInput(:,t==0,:,:);
%theoreticalInput(:,t<0,:,:) = 0;
if centerData.On
    theoreticalInput(1,:,:,:) = theoreticalInput(1,:,:,:) - theoreticalInput(1,:,centerData.inds(1),centerData.inds(2));
    theoreticalInput(2,:,:,:) = theoreticalInput(2,:,:,:) - theoreticalInput(2,:,centerData.inds(1),centerData.inds(2));
else
    theoreticalInput(1,:,:,:) = theoreticalInput(1,:,:,:);% - theoreticalInput(1,:,2,2);
    theoreticalInput(2,:,:,:) = theoreticalInput(2,:,:,:);% - theoreticalInput(2,:,3,3);    
end
theoreticalInput = theoreticalInput./max(abs(theoreticalInput),[],[2,3,4]);

%% Fit dynamic model
switch fitType
    case 'sigmasHeldConstant'
        N = size(RtoFit,4);
        M = size(theoreticalInput,1);
        
        overlapsSet = zeros(N,N+M);
        overlapsSet(1,1) = 1;
        overlapsSet(2,1) = -1.7;
        overlapsSet(2,2) = 0.5;
        overlapsSet(1,N+1) = 0.5;
        overlapsSet(1,N+2) = 0;
        overlapsSet(2,N+2) = 0;
        overlapsSet(2,N+1) = 1.9;
        
        % P0 = reshape(overlapsSet,[N*(N+M),1]);
        if any(isnan(P0))
            P0 = zeros(N*(N+M),1);
        end            
        lb = -Inf*ones(N,N+M);
        temp = zeros(N);
        temp(~eye(N)) = -Inf*1;
        lb(1:N,1:N) = temp;
        lb = reshape(lb,[N*(N+M),1]);
        % lb = -Inf*ones(N*(N+M),1);
        % lb(3) = -Inf;
        % lb(N*(N+M)-1) = -Inf;
        % lb(N*(N+M)) = -Inf;
        % ub = Inf*ones(N*(N+M),1);
        % ub(3) = Inf;
        % ub(N*(N+M)-1) = Inf;
        % ub(N*(N+M)) = Inf;
        
        if any(isnan(P0))
            P0 = zeros(N*(N+M),1);
        end
        if any(isnan(lb))
            lb = -Inf*ones(N*(N+M),1);
        end
        if any(isnan(ub))
            ub = Inf*ones(N*(N+M),1);
        end
        sigmas = ones(N+M,1);
        
        OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10,'Display','iter',...
            'StepTolerance',1e-40);
        R0 = nanmean(RtoFit(1,:,:,:),[2,3]);
        minimizer = @(P)minimizant(P,R0,theoreticalInput,RtoFit,t,tau,sigmas);
        [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
        
        modelFEF.tau = tau;
        modelFEF.R0 = R0;
        modelFEF.overlaps = reshape(P,[N,N+M]);
        modelFEF.sigmas = sigmas;
        modelFEF.fval = fval;
        modelFEF.exitflag = exitflag;
        modelFEF.output = output;
        modelFEF.lambda = lambda;
        modelFEF.grad = grad;
        modelFEF.hessian = hessian;
        modelFEF.t = t;
        modelFEF.dt = dt;
        
    case 'fitSigmas'
        
        % Fit sigmas
        N = size(RtoFit,4);
        M = size(theoreticalInput,1);
        
        overlapsSet = zeros(N,N+M);
        overlapsSet(1,1) = 1;
        overlapsSet(2,1) = -1.7;
        overlapsSet(2,2) = 0.5;
        overlapsSet(1,N+1) = 0.5;
        overlapsSet(1,N+2) = 0;
        overlapsSet(2,N+2) = 0;
        overlapsSet(2,N+1) = 1.9;
        
        % P0 = reshape(overlapsSet,[N*(N+M),1]);
        if any(isnan(P0))
            P0 = zeros(N*(N+M) + N + M,1);    
        end
        lb = -Inf*ones(N,N+M);
        temp = zeros(N);
        temp(~eye(N)) = -Inf*1;
        lb(1:N,1:N) = temp;
        lb = reshape(lb,[N*(N+M),1]);
        lb = [lb; zeros(N+M,1)];
        ub = Inf*ones(N*(N+M)+N+M,1);
        
        %lb(end-(N+M)-1) = 0;
        %ub(end-(N+M)-1) = 0;
        %lb(end-(N+M)) = 0;
        %ub(end-(N+M)) = 0;
        % lb = -Inf*ones(N*(N+M),1);
        % lb(3) = -Inf;
        % lb(N*(N+M)-1) = -Inf;
        % lb(N*(N+M)) = -Inf;
        % ub = Inf*ones(N*(N+M),1);
        % ub(3) = Inf;
        % ub(N*(N+M)-1) = Inf;
        % ub(N*(N+M)) = Inf;
        
        if any(isnan(P0))
            P0 = zeros(N*(N+M),1);
        end
        if any(isnan(lb))
            lb = -Inf*ones(N*(N+M),1);
        end
        if any(isnan(ub))
            ub = Inf*ones(N*(N+M)+N+M,1);
        end
        
        OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10,'Display','iter',...
            'StepTolerance',1e-40);
        R0 = nanmean(RtoFit(1,:,:,:),[2,3]);
        minimizer = @(P)minimizantSigmas(P,R0,theoreticalInput,RtoFit,t,tau);
        [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
        
        modelFEF.tau = tau;
        modelFEF.R0 = R0;
        modelFEF.overlaps = reshape(P(1:end-(N+M)),[N,N+M]);
        modelFEF.sigmas = P(end-(N+M)+1:end);
        modelFEF.fval = fval;
        modelFEF.exitflag = exitflag;
        modelFEF.output = output;
        modelFEF.lambda = lambda;
        modelFEF.grad = grad;
        modelFEF.hessian = hessian;
        modelFEF.t = t;
        modelFEF.dt = dt;
        modelFEF.dimNames = dimNames;

    case 'fitSigmasWithDummy'
        
        % Fit sigmas
        N = size(RtoFit,4)+1;
        M = size(theoreticalInput,1);
        
        if isempty(initialConditions)
            P0 = zeros(N*(N+M) + N + M,1);            
        else
            P0 = initialConditions;
        end
        
        lb = -Inf*ones(N,N+M);
        temp = zeros(N);
        temp(~eye(N)) = -Inf*1;
        lb(1:N,1:N) = temp;
        lb(N,:) = 0;
        lb = reshape(lb,[N*(N+M),1]);
        lb = [lb; zeros(N+M,1)];
        
        ub = Inf*ones(N,N+M);
        ub(N,:) = 0;
        ub(N,N) = Inf;
        ub(N,N+M) = Inf;
        ub = reshape(ub,[N*(N+M),1]);
        ub = [ub; Inf*ones(N+M,1)];
                        
        if any(isnan(P0))
            P0 = zeros(N*(N+M),1);
        end
        if any(isnan(lb))
            lb = -Inf*ones(N*(N+M),1);
        end
        if any(isnan(ub))
            ub = Inf*ones(N*(N+M)+N+M,1);
        end
        
        OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10,'Display','iter',...
            'StepTolerance',1e-40);
        R0 = nanmean(RtoFit(1,:,:,:),[2,3]);
        minimizer = @(P)minimizantSigmasDummy(P,R0,theoreticalInput,RtoFit,t,tau,1);
        [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
        
        modelFEF.tau = tau;
        modelFEF.R0 = [squeeze(R0); 0];
        modelFEF.overlaps = reshape(P(1:end-(N+M)),[N,N+M]);
        modelFEF.sigmas = P(end-(N+M)+1:end);
        modelFEF.fval = fval;
        modelFEF.exitflag = exitflag;
        modelFEF.output = output;
        modelFEF.lambda = lambda;
        modelFEF.grad = grad;
        modelFEF.hessian = hessian;
        modelFEF.t = t;
        modelFEF.dt = dt;
        modelFEF.dimNames = dimNames;
        
    case 'fitSigmasExtraInput'
        
        % Fit sigmas
        theoreticalInput(size(theoreticalInput,1)+1,:,:,:) = 0;
        N = size(RtoFit,4);
        M = size(theoreticalInput,1);
        
        
        if any(isnan(P0))
            P0 = zeros(N*(N+M) + N + M,1);
            P0 = [P0; 1/100; 250];
        end
        lb = -Inf*ones(N,N+M);
        temp = zeros(N);
        temp(~eye(N)) = -Inf*1;
        lb(1:N,1:N) = temp;
        lb = reshape(lb,[N*(N+M),1]);
%         lb = [lb; zeros(N+M,1); 0; -Inf];
        lb = [lb; sigmaLowerBound*ones(N+M,1); 0; -Inf];
        ub = Inf*ones(N*(N+M)+N+M,1);
        ub = [ub; 1/200; Inf];
        
        
        if any(isnan(P0))
            P0 = zeros(N*(N+M),1);
        end
        if any(isnan(lb))
            lb = -Inf*ones(N*(N+M),1);
        end
        if any(isnan(ub))
            ub = Inf*ones(N*(N+M)+N+M,1);
        end
        
        OPTIONS = optimoptions('fmincon','MaxFunEvals',3e10,'MaxIterations',1e10,'Display','iter',...
            'StepTolerance',1e-40);
        R0 = nanmean(RtoFit(1,:,:,:),[2,3]);
        minimizer = @(P)minimizantSigmasExtraInput(P,R0,theoreticalInput,RtoFit,t,tau,3);
        [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
        
        Ptemp = P;
        theoreticalInput(3,:,:,:) = repmat(1./(1 + exp( -(t-Ptemp(end))*Ptemp(end-1) )),[1,1,size(theoreticalInput,3),size(theoreticalInput,4)]);
        Ptemp = P(1:end-2);
        
        modelFEF.tau = tau;
        modelFEF.R0 = R0;
        modelFEF.overlaps = reshape(Ptemp(1:end-(N+M)),[N,N+M]);
        modelFEF.sigmas = Ptemp(end-(N+M)+1:end);
        modelFEF.fval = fval;
        modelFEF.exitflag = exitflag;
        modelFEF.output = output;
        modelFEF.lambda = lambda;
        modelFEF.grad = grad;
        modelFEF.hessian = hessian;
        modelFEF.t = t;
        modelFEF.dt = dt;
        modelFEF.dimNames = dimNames;
        modelFEF.extraInput = P(end-1:end);
        
    otherwise
        error(['fitType ' fitType ' not recognized!'])
end

%% probe the objective function by perturbations around the found minimum
if probeObjectiveFunction
    [D12,D21] = meshgrid(linspace(-2,2,11),linspace(-2,2,11));
    del = zeros([size(modelFEF.overlaps),numel(D12)]);
%     del(2,3,:) = -1;
    del(1,2,:) = D12(:);
    del(2,1,:) = D21(:);
    err = probeObjective(minimizer,modelFEF,del);
    figure
    surf(D12,D21,reshape(log(err),size(D12)),'EdgeColor','none');
    xlabel('x')
end

%% Make predictions
overlapsTemp = modelFEF.overlaps;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',theoreticalInput(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

Rhat = permute(kappas,[2,3,4,1]);
%% Compare fit quality
for si = 1:length(speedsFEF)
    rvalAll = corrcoef([squeeze(R_n(taccept,si,2,:)) squeeze(Rhat(taccept,si,2,:))]);
    rvals(:,si) = diag(rvalAll(size(R,4)+1:end,1:size(R,4)));
end
rvalsZ = 0.5*log( (1+rvals)./(1-rvals) );
[~,rvalsort] = sort(mean(rvals,2));

RhatTemp = Rhat(taccept,:,:,:);
Yhat = nan(size(RtoFit,1)*numel(trainCondition(:)),size(RtoFit,4));
Y = nan(size(RtoFit,1)*numel(trainCondition(:)),size(RtoFit,4));
conditions = nan(size(RtoFit,1)*numel(trainCondition(:)),3);

idx = 1;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        Y(idx:idx+size(RtoFit,1)-1,:) = squeeze(RtoFit(:,si,ci,:));
        Yhat(idx:idx+size(RtoFit,1)-1,:) = squeeze(Rhat(:,si,ci,1:size(RtoFit,4)));
        conditions(idx:idx+size(RtoFit,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit,1),1]);
        idx = idx + size(RtoFit,1);
    end
end
rvalsFit = corrcoef([Y, Yhat]);
RhatCC = diag(rvalsFit(size(Y,2)+1:end,1:size(Y,2)));
RhatZ = 0.5*log( (1+RhatCC)./(1-RhatCC) );
[~,zvalSort] = sort(RhatZ);
    
end

%% Saving
if saveResults
        
    if isempty(saveLocation)
        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' subject ...
            '/ReducedRankModel'];
        if ~exist(saveLocation,'dir')
            mkdir(saveLocation)
        end
        if reducedSave
            save([saveLocation '/fitReducedRankModel' datestr(now,'yyyymmdd')],'-v7.3',...
                'N','speeds','meanEyeSpeed','eye_t','Rhat','t','initGain','inputTheoretical','modelFEF','RtoFit')
        else
            save([saveLocation '/fitReducedRankModel' datestr(now,'yyyymmdd')],'-v7.3')
        end
    else
        if reducedSave
            save(saveLocation,'-v7.3',...
                'N','speeds','meanEyeSpeed','eye_t','Rhat','t','initGain','inputTheoretical','modelFEF','RtoFit')
        else
            save(saveLocation,'-v7.3')
        end
    end
    
end

%% Plotting
if plotOpts.On

    %% Fit quality
    figure
    for neuroni = 1:size(RtoFit,4)
        tempFEF = RtoFit(:,:,:,neuroni);
        test = Rhat(:,:,:,neuroni);
        for ci = 1:length(cohsFEF)
            for si = 1:length(speedsFEF)
                subplot(length(cohsFEF),length(speedsFEF),si + (ci-1)*length(speedsFEF))
                plot(tempFEF(taccept,si,ci),test(taccept,si,ci),'o','Color',speedColors(si,:),...
                    'DisplayName',['Dim = ' num2str(neuroni) ', Speed = ' num2str(speedsFEF(si)) ', Coh = ' num2str(cohsFEF(ci))])
                hold on
            end
        end
    end
    
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(cohsFEF),length(speedsFEF),si + (ci-1)*length(speedsFEF))
            axis square
            plotUnity;
            xlabel('FEF data')
            ylabel('Fit data')
        end
    end
    
    %% Fit quality distribution
    figure('Name','Fit quality distribution')
    plot(RhatZ,'bo')
    hold on
    xlabel('Dimension')
    ylabel('Fisher transformed r-values')
        
    %% Model parameters
    hModelParameters = figure('Name','ModelParameters');
    subplot(2,2,1)
    imagesc(modelFEF.overlaps(:,1:N))
    colorbar
    axis square
    ylabel('Reduced rank dimension')
    xlabel('Reduced rank dimension')
    
    subplot(2,2,2)
    imagesc(modelFEF.overlaps(:,N+1:end))
    colorbar
    axis equal
    axis tight
    ylabel('Reduced rank dimension')
    xlabel('Input dimension')
    
    subplot(2,2,3)
    temp = nan(N^2,1);
    tickNames = {};
    for i = 1:N
        temp(i) = modelFEF.overlaps(i,i);
%         tickNames{i} = [sprintf(strrep('\u03C3', '\u', '\x')) '_{' modelFEF.dimNames{i} ',' modelFEF.dimNames{i} '}'];
        if i <= length(modelFEF.dimNames)
            tickNames{i} = [modelFEF.dimNames{i} ',' modelFEF.dimNames{i}];
        else
            tickNames{i} = '';
        end
    end
    ind = i;
    for i = 1:N
        for j = 1:N
            if i ~= j
                ind = ind+1;
                temp(ind) = modelFEF.overlaps(i,j);
%                 tickNames{ind} = [sprintf(strrep('\u03C3', '\u', '\x')) '_{' modelFEF.dimNames{i} ',' modelFEF.dimNames{j} '}'];
                if i <= length(modelFEF.dimNames) && j <= length(modelFEF.dimNames)
                    tickNames{ind} = [modelFEF.dimNames{i} ',' modelFEF.dimNames{j}];
                else
                    tickNames{ind} = '';
                end
            end
        end
    end
    bar(temp,'DisplayName',[subject ' ' fitType ' Recurrent connections'])
    axTemp = gca;
    axTemp.XTickLabels = tickNames;
    axTemp.YLabel.String = 'Overlap';
    %text(min(xlim),max(ylim),tickNames)
    
    
    subplot(2,2,4)
    temp = nan(M*N,1);
    tickNames = {};
    ind = 0;
    for i = 1:N
        for j = 1:M
            ind = ind+1;
            temp(ind) = modelFEF.overlaps(i,N+j);
%             tickNames{ind} = [sprintf(strrep('\u03C3', '\u', '\x')) '_{' modelFEF.dimNames{i} ', Input' num2str(j) '}'];
            if i <= length(modelFEF.dimNames)
                tickNames{ind} = [modelFEF.dimNames{i} ', Input' num2str(j)];
            else
                tickNames{ind} = [' , Input' num2str(j)];
            end
        end
    end
    bar(temp,'DisplayName',[subject ' ' fitType ' Recurrent connections'])
    axTemp = gca;
    axTemp.XTickLabels = tickNames;
    axTemp.YLabel.String = 'Overlap';
    %text(min(xlim),max(ylim),tickNames)
    
    %% 
            
        hModelandData = figure('Name','Predicted response along targeted dimensions','Position',[63 169 1606 1079]);
        for dimi = 1:M
            subplot(N+M,2,2*(dimi-1)+1)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,theoreticalInput(dimi,:,speedi,cohi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100],...
                        'DisplayName',['Speed = ' num2str(speedsFEF(speedi)) ', Coh = ' num2str(cohsFEF(cohi))])
                    hold on
                end
            end
            axis tight
            title('Input from MT')
            xlabel('Time from motion onset (ms)')
            ylabel(['Input ' num2str(dimi)])
        end
            
        
        for dimi = 1:N
            subplot(N+M,2,2*(dimi-1)+2*M+1)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,Rhat(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100],...
                        'DisplayName',['Speed = ' num2str(speedsFEF(speedi)) ', Coh = ' num2str(cohsFEF(cohi))])
                    hold on
                end
            end
            axis tight
            title('Fit predictions')
            xlabel('Time from motion onset (ms)')
            if dimi > length(dimNames)
                ylabel(['Dummy related activity (a.u.)'])
            else
                ylabel([dimNames{dimi} ' related activity (a.u.)'])
            end
                        
            subplot(N+M,2,2*(dimi-1)+2*M+2)
            if dimi <= length(dimNames)
                for speedi = 1:length(speedsFEF)
                    for cohi = 1:length(cohsFEF)
                        plot(t,RtoFit(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100],...
                            'DisplayName',['Speed = ' num2str(speedsFEF(speedi)) ', Coh = ' num2str(cohsFEF(cohi))])
                        hold on
                    end
                end
                axis tight
                title('Data')
                xlabel('Time from motion onset (ms)')
                ylabel([dimNames{dimi} ' related activity (a.u.)'])
            end
        end
        
        compLims = repmat([Inf,-Inf],[N,1]);
        for dimi = 1:N
            for compi = 1:2
                subplot(N+M,2,2*(dimi-1)+2*M+compi)
                limTemp = ylim;
                compLims(dimi,1) = min([compLims(dimi,1),limTemp(1)]);
                compLims(dimi,2) = max([compLims(dimi,2),limTemp(2)]);
            end
        end
        for dimi = 1:N
            for compi = 1:2
                subplot(N+M,2,2*(dimi-1)+2*M+compi)
                ylim(compLims(dimi,:))
            end
        end
        
        figure('Name','Systematic differencs from observed targeted dimensions','Position',[63 169 1606 1079])
        for dimi = 1:2
            subplot(3,1,dimi)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,RtoFit(:,speedi,cohi,dimi)-Rhat(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            title('Difference from fit')
            xlabel('Time from motion onset (ms)')
            ylabel([dimNames{dimi} ' related activity (a.u.)'])
        end
        
        gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1456 59 1070 1263]);
        for targetedDim = 1:N
            tempData = [];
            subplot(N,2,2*targetedDim-1)            
            for si = 1:length(speeds)
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
                initRatesTemp = squeeze(nanmean(Rhat(t >= 700 & t <= 800,si,:,targetedDim),1));
                plot(eyeSpeedTemp(si,:),initRatesTemp,...
                    'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                    'DisplayName',['Speed = ' num2str(speeds(si))])
                hold on
                tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
            end
            for si = 1:length(speeds)
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
                initRatesTemp = squeeze(nanmean(Rhat(t >= 100 & t <= 200,si,:,targetedDim),1));
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
            if targetedDim <= length(dimNames)
                ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
            end
            
            tempData = [];
            subplot(N,2,2*targetedDim)
            for si = 1:length(speedsFEF)
                initRatesTemp = squeeze(nanmean(Rhat(t >= 700 & t <= 800,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                    'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                    'DisplayName',['Steady State pursuit, Speed = ' num2str(speedsFEF(speedi))])
                hold on
                tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
                initRatesTemp = squeeze(nanmean(Rhat(t >= 100 & t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                    'DisplayName',['Pursuit initiation, Speed = ' num2str(speedsFEF(speedi))])
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
            gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
            axis square
            ax = axis;
            text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
            xlabel('Behavioral gain (unitless)')
            if targetedDim <= length(dimNames)
                ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
            end
        end
    
end

%% Save figures
if saveFigures
    
    if isempty(saveLocation)
        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' subject ...
            '/fitReducedRankModel/' datestr(now,'yyyymmdd')];
    end
        
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    
    savefig(hModelandData,[saveLocation '/reducedRankModelVsTargetedDimensions.fig'])
    savefig(gvth ,[saveLocation '/reducedRankModelvsGain.fig'])
    savefig(hModelParameters ,[saveLocation '/reducedRankModelParameters.fig'])
    
    
end

%% Functions
    
 %% Objective
function out = minimizant(P,R0,inputs,R,t,tau,sigmas)
    N = size(R,4);
    M = size(inputs,1);
    overlaps = reshape(P,[N N+M]);
    
    for si = 1:size(inputs,3)
        for ci = 1:size(inputs,4)
            [~,kappas(:,:,si,ci),sigmaTildes(:,:,:,si,ci),deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/(t(2)-t(1)),...
                't',t,...
                'us',inputs(:,:,si,ci),...
                'kappas0',squeeze(R0),...
                'overlaps',overlaps,...
                'sigmas',sigmas);
            
        end
    end
    Rest = permute(kappas,[2,3,4,1]);
        
    out = nansum((R(:) - Rest(:)).^2);
    
%% Include sigmas in fit
function out = minimizantSigmas(P,R0,inputs,R,t,tau)
    N = size(R,4);
    M = size(inputs,1);
    sigmas = P(end-(N+M)+1:end);
    P = P(1:end-(N+M));
    overlaps = reshape(P,[N N+M]);
    
    for si = 1:size(inputs,3)
        for ci = 1:size(inputs,4)
            [~,kappas(:,:,si,ci),sigmaTildes(:,:,:,si,ci),deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/(t(2)-t(1)),...
                't',t,...
                'us',inputs(:,:,si,ci),...
                'kappas0',squeeze(R0),...
                'overlaps',overlaps,...
                'sigmas',sigmas);
            
        end
    end
    Rest = permute(kappas,[2,3,4,1]);
        
    out = nansum((R(:) - Rest(:)).^2);
    
%% Include sigmas in fit w/ dummy population mode
function out = minimizantSigmasDummy(P,R0,inputs,R,t,tau,dummyModes)
    N = size(R,4)+dummyModes;
    M = size(inputs,1);
    sigmas = P(end-(N+M)+1:end);
    P = P(1:end-(N+M));
    overlaps = reshape(P,[N N+M]);
    kappas0 = [squeeze(R0); zeros(dummyModes,1)];
    
    for si = 1:size(inputs,3)
        for ci = 1:size(inputs,4)
            [~,kappas(:,:,si,ci),sigmaTildes(:,:,:,si,ci),deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/(t(2)-t(1)),...
                't',t,...
                'us',inputs(:,:,si,ci),...
                'kappas0',kappas0,...
                'overlaps',overlaps,...
                'sigmas',sigmas);
            
        end
    end
    Rest = permute(kappas(1:size(R,4),:,:,:),[2,3,4,1]);
        
    out = nansum((R(:) - Rest(:)).^2);
    
%% Include sigmas and parameters for additional input in fit 
function out = minimizantSigmasExtraInput(P,R0,inputs,R,t,tau,extraInd)
    N = size(R,4);
    M = size(inputs,1);
    Ptemp = P;
    inputParams = Ptemp(end-1:end);
    Ptemp = Ptemp(1:end-2);
    sigmas = Ptemp(end-(N+M)+1:end);
    Ptemp = Ptemp(1:end-(N+M));
    overlaps = reshape(Ptemp,[N N+M]);
    kappas0 = squeeze(R0);
    
    inputs(extraInd,:,:,:) = repmat(1./(1 + exp( -(t-inputParams(2))*inputParams(1) )),[1, 1, size(inputs,3), size(inputs,4)]);
    
    for si = 1:size(inputs,3)
        for ci = 1:size(inputs,4)
            [~,kappas(:,:,si,ci),sigmaTildes(:,:,:,si,ci),deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/(t(2)-t(1)),...
                't',t,...
                'us',inputs(:,:,si,ci),...
                'kappas0',kappas0,...
                'overlaps',overlaps,...
                'sigmas',sigmas);
            
        end
    end
    Rest = permute(kappas(1:size(R,4),:,:,:),[2,3,4,1]);
        
    out = nansum((R(:) - Rest(:)).^2);

%% outputsObjective
function out = minimizantOutputs(P,R0,inputs,R,t,tau,sigmas,initGain,speedsFEF)
    N = size(R,4);
    M = size(inputs,1);
    overlaps = reshape(P,[N N+M]);
    
    for si = 1:size(inputs,3)
        for ci = 1:size(inputs,4)
            [~,kappas(:,:,si,ci),sigmaTildes(:,:,:,si,ci),deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/(t(2)-t(1)),...
                't',t,...
                'us',inputs(:,:,si,ci),...
                'kappas0',squeeze(R0),...
                'overlaps',overlaps,...
                'sigmas',sigmas);
            
        end
    end
    Rest = permute(kappas,[2,3,4,1]);
    
    tempGains = [];
    tempSpeeds = [];
    for si = 1:length(speedsFEF)
        initRatesTemp = squeeze(nanmean(Rest(t >= 700 & t <= 800,si,:,2),1));
        tempGains = [tempGains; initGain(si,:,3)',initRatesTemp(:)];
        initRatesTemp = squeeze(nanmean(Rest(t >= 100 & t <= 200,si,:,2),1));
        tempGains = [tempGains; initGain(si,:,2)',initRatesTemp(:)];
        
        initRatesTemp = squeeze(nanmean(Rest(t >= 700 & t <= 800,si,:,1),1));
        tempSpeeds = [tempGains; speedsFEF(:), initRatesTemp(:)];
        initRatesTemp = squeeze(nanmean(Rest(t >= 100 & t <= 200,si,:,1),1));
        tempSpeeds = [tempGains; speedsFEF(:), initRatesTemp(:)];
    end
    cc_gain = corrcoef(tempGains(:,1),tempGains(:,2));
    cc_speed = corrcoef(tempSpeeds(:,1),tempSpeeds(:,2));
        
    out = -cc_speed - cc_gain;
    
%% probeObjective
function err = probeObjective(minimizer,modelFEF,del)
%%
    N = size(modelFEF.overlaps,1);
    M = size(modelFEF.overlaps,2)-N;
    for ind = 1:size(del,3)

        overlapsTemp = modelFEF.overlaps;
        overlapsTemp = overlapsTemp + del(:,:,ind);

        err(ind) = minimizer([reshape(overlapsTemp,[N*(N+M),1]); modelFEF.sigmas]);

    end
