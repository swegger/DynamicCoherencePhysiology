function reducedRankDynamicsModel(varargin)
%% MTtoFEFregression
%
%
%%

%% Defaults
plotOpts_default.On = true;
speedPrefOpts_default.tWin = [40,120];
speedPrefOpts_default.P0 = [16,1];
speedPrefOpts_default.ub = [128,128];
speedPrefOpts_default.lb = [0, 0];
speedPrefOpts_default.c = NaN;
speedPrefOpts_default.s = NaN;
speedPrefOpts_default.d = 0;
theoretical_default.weightTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';
simulateMT_default.On = false;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'tau',20)
addParameter(Parser,'kappas0',NaN)
addParameter(Parser,'overlaps',NaN)
addParameter(Parser,'sigmas',NaN)
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'simulateMT',simulateMT_default)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,varargin{:})

tau = Parser.Results.tau;
kappas0 = Parser.Results.kappas0;
overlaps = Parser.Results.overlaps;
sigmas = Parser.Results.sigmas;
objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
simulateMT = Parser.Results.simulateMT;
theoretical = Parser.Results.theoretical;
plotOpts = Parser.Results.plotOpts;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
if simulateMT.On
    t = -500:900;
    speeds = speedsFEF;
    [MT, spref, cohs] = simulateMTdata(mtResults,mtNeuron_t,simulateMT.modelN,...
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
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);
inputs_z = (inputs-mean(inputs,[2,3,4]))./...
    std(inputs,[],[2,3,4]);
inputs_n = inputs./max(abs(inputs),[],[2,3,4]);

sp = spref(~isnan(interpolatedR(1,1,1,:)));
sw = swidth(~isnan(interpolatedR(1,1,1,:)));
%% A theoretical treatment, rather than data driven

switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(sp)'-log2(mean(speedsFEF))) ones(size(sp'))];
        
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

% Normalize 
Atheory = Atheory./vecnorm(Atheory);

% Compute output along theoretical channels
for ri = 1:size(Atheory,2)
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            theoreticalInput(ri,:,si,ci) = inputs_n(:,:,si,ci)'*Atheory(:,ri);
        end
    end
end
theoreticalInput = theoreticalInput - theoreticalInput(:,t==0,:,:);
theoreticalInput(:,t<0,:,:) = 0;
theoreticalInput(1,:,:,:) = theoreticalInput(1,:,:,:) - theoreticalInput(1,:,2,2);
theoreticalInput(2,:,:,:) = theoreticalInput(2,:,:,:) - theoreticalInput(2,:,3,3);
theoreticalInput = theoreticalInput./max(abs(theoreticalInput),[],[2,3,4]);
theoreticalInput(3,:,:,:) = 0;
% theoreticalInput(2,:,:,:) = theoreticalInput(2,:,:,:)/10000;

%% Set overlaps and sigmas
N = 2;
M = size(theoreticalInput,1);
if any(isnan(overlaps(:)))
    overlapsSet = zeros(N,N+M);
    overlapsSet(1,1) = 1;
    overlapsSet(2,1) = 1.7;
    overlapsSet(2,2) = 0.5;
    overlapsSet(1,N+1) = 0.5;
    overlapsSet(1,N+2) = 0;
    overlapsSet(2,N+2) = 0;
    overlapsSet(2,N+1) = -1.9;
else
    overlapsSet = overlaps;
end

if any(isnan(sigmas))
    sigmasSet = ones(N+M,1)*1; %randn(N+M,1);
else
    sigmasSet = sigmas;
end

%% Simulate reduced rank network
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci),sigmaTildes(:,:,:,si,ci),deltas(:,si,ci),~] = simulateLatentDynamics('tau',tau/(t(2)-t(1)),...
            't',t,...
            'us',theoreticalInput(:,:,si,ci),...
            'kappas0',kappas0,...
            'overlaps',overlapsSet,...
            'sigmas',sigmasSet);
            
    end
end

% Calculate gain
g = sigmaTildes./overlapsSet;

%% Plotting
if plotOpts.On
    
    %% Plot reduced rank outputs
    dims = {'speed','gain'};
    figure('Name',['Inputs and outputs'],'Position',[1961 282 560 1040])
    for di = 1:M
        subplot(length(dims)+M,2,1+(di-1)*2)
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                plot(t,theoreticalInput(di,:,si,ci),...
                    'Color',[speedColors(si,:) cohsFEF(ci)/100])
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Input ' num2str(di)])
    end
    for di = 1:length(dims)
        subplot(length(dims)+M,2,2*M+(di-1)*2+1)
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                plot(t,kappas(di,:,si,ci),...
                    'Color',[speedColors(si,:) cohsFEF(ci)/100])
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['\kappa_{' dims{di} '}'])
    end
    
    for di = 1:length(dims)
        subplot(length(dims)+M,2,2*M+2+(di-1)*2)
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                plot(t,squeeze(deltas(:,si,ci)),...
                    'Color',[speedColors(si,:) cohsFEF(ci)/100])
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['\Delta'])
    end
    
end