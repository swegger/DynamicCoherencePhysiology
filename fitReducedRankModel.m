function modelFEF = fitReducedRankModel(nerualDataFile,varargin)
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
equalizeInputsPriorToStimulusOnset_default.On = false;
simulateMT_default.On = false;

centerData_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'nerualDataFile')
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
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,nerualDataFile,varargin{:})

nerualDataFile = Parser.Results.nerualDataFile;
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
saveOpts = Parser.Results.saveOpts;

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
temp = load(nerualDataFile,'BinitOrth','pInit','initCoh','initGain');

BinitOrth = temp.BinitOrth;
pInit = temp.pInit;
fef_t = temp.initCoh.neuron_t;
initGain = temp.initGain;

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
if simulateMT.On
    t = -500:900;
    speeds = speedsFEF;
    [MT, spref, cohs] = simulateMTdata(mtResults,t,simulateMT.modelN,...
        'speedsMT',speeds,'removeBaseline',simulateMT.removeBaseline,'gaussianApprox',simulateMT.gaussianApprox);
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
R = pInit(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);
dt = t(2)-t(1);
t = t(iMT);

inputs_z = (inputs-mean(inputs(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(inputs(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);
inputs_n = inputs./max(abs(inputs),[],[2,3,4]);

R_z = (R-mean(R,[1,2,3]))./std(R,[],[1,2,3]);
R_n = R./max(abs(R),[],[1,2,3]);

taccept = t >= tWin(1) & t <= tWin(2);

RtoFit = R_n(:,:,:,[1,2,3]);
RtoFit = RtoFit - mean(RtoFit(t<=0,:,:,:),[1,2,3]);
if centerData.On
    RtoFit = RtoFit - RtoFit(:,centerData.inds(1),centerData.inds(2),:);
    initGain = initGain - initGain(centerData.inds(1),centerData.inds(2),:);
end
RtoFit = RtoFit./max(abs(RtoFit),[],[1,2,3]);

%% Set MT weights
switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
%         Atheory = [(log2(sp)' - log2(mean(speedsFEF))) ones(size(sp'))];
        Atheory = [(log2(sp)') ones(size(sp'))];
        
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
theoreticalInput(3,:,:,:) = repmat(permute(mean(R_n(:,:,:,4),[2,3]),[2,1]),[1,1,length(speedsFEF),length(cohsFEF)]);
theoreticalInput = theoreticalInput - theoreticalInput(:,t==0,:,:);
%theoreticalInput(:,t<0,:,:) = 0;
if centerData.On
    theoreticalInput(1,:,:,:) = theoreticalInput(1,:,:,:) - theoreticalInput(1,:,centerData.inds(1),centerData.inds(2));
    theoreticalInput(2,:,:,:) = theoreticalInput(2,:,:,:) - theoreticalInput(2,:,centerData.inds(1),centerData.inds(2));
else
    theoreticalInput(1,:,:,:) = theoreticalInput(1,:,:,:) - theoreticalInput(1,:,2,2);
    theoreticalInput(2,:,:,:) = theoreticalInput(2,:,:,:) - theoreticalInput(2,:,3,3);    
end
theoreticalInput = theoreticalInput./max(abs(theoreticalInput),[],[2,3,4]);

%% Fit dynamic model
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
P0 = zeros(N*(N+M),1);
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

%% probe the objective function by perturbations around the found minimum
if probeObjectiveFunction
    [D14,D31] = meshgrid(linspace(-2,2,11),linspace(0,4,11));
    del = zeros([size(modelFEF.overlaps),numel(D14)]);
%     del(2,3,:) = -1;
    del(1,4,:) = D14(:);
    del(3,1,:) = D31(:);
    err = probeObjective(minimizer,modelFEF,del);
    figure
    surf(D14,D31,reshape(err,size(D14)),'EdgeColor','none');
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
        Yhat(idx:idx+size(RtoFit,1)-1,:) = squeeze(Rhat(:,si,ci,:));
        conditions(idx:idx+size(RtoFit,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit,1),1]);
        idx = idx + size(RtoFit,1);
    end
end
rvalsFit = corrcoef([Y, Yhat]);
RhatCC = diag(rvalsFit(size(Y,2)+1:end,1:size(Y,2)));
RhatZ = 0.5*log( (1+RhatCC)./(1-RhatCC) );
[~,zvalSort] = sort(RhatZ);

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
                plot(tempFEF(taccept,si,ci)*1000,test(taccept,si,ci)*1000,'o','Color',speedColors(si,:))
                hold on
            end
        end
    end
    
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(cohsFEF),length(speedsFEF),si + (ci-1)*length(speedsFEF))
            axis square
            plotUnity;
            xlabel('FEF data (sp/s)')
            ylabel('Fit data (sp/s)')
        end
    end
    
    %% Fit quality distribution
    figure('Name','Fit quality distribution')
    plot(RhatZ(zvalSort),'bo')
    hold on
    xlabel('FEF neuron #')
    ylabel('Fisher transformed r-values')
        
    %% Model parameters
    figure
    subplot(1,2,1)
    imagesc(modelFEF.overlaps(:,1:N))
    colorbar
    axis square
    ylabel('Reduced rank dimension')
    xlabel('Reduced rank dimension')
    
    subplot(1,2,2)
    imagesc(modelFEF.overlaps(:,N+1:end))
    colorbar
    axis equal
    axis tight
    ylabel('Reduced rank dimension')
    xlabel('Input dimension')
    
    %% 
    
    
        dimNames = {'Speed','Coherence','Gain','Offset'};
        
        figure('Name','Predicted response along targeted dimensions','Position',[63 169 1606 1079])
        for dimi = 1:M
            subplot(N+M,2,2*(dimi-1)+1)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,theoreticalInput(dimi,:,speedi,cohi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            xlabel('Time from motion onset (ms)')
            ylabel(['Input ' num2str(dimi)])
        end
            
        
        for dimi = 1:N
            subplot(N+M,2,2*(dimi-1)+2*M+1)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,Rhat(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            title('Fit predictions')
            xlabel('Time from motion onset (ms)')
            ylabel([dimNames{dimi} ' related activity (a.u.)'])
                        
            subplot(N+M,2,2*(dimi-1)+2*M+2)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(t,RtoFit(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                    hold on
                end
            end
            axis tight
            title('Data')
            xlabel('Time from motion onset (ms)')
            ylabel([dimNames{dimi} ' related activity (a.u.)'])
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
        
        gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
        for targetedDim = 1:N
            tempData = [];
            subplot(N,1,targetedDim)
            for si = 1:length(speedsFEF)
                initRatesTemp = squeeze(nanmean(Rhat(t >= 700 & t <= 800,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                    'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                hold on
                tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
                initRatesTemp = squeeze(nanmean(Rhat(t >= 100 & t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
            gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
            axis square
            ax = axis;
            text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
            xlabel('Behavioral gain (unitless)')
            ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
        end
    
end


%% Saving
if saveOpts.On
    save(saveOpts.location);
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

        err(ind) = minimizer(reshape(overlapsTemp,[N*(N+M),1]));

    end
