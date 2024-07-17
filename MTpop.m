function [inputs, inputs_n, theoreticalInput, theoreticalInputNull, spref] = MTpop(varargin)
%%
%
%
%
%
%%

%% Defaults
simulateMT_default.On = false;
theoretical_default.weightTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';
equalizeInputsPriorToStimulusOnset_default.On = false;
centerData_default.On = false;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'simulateMT',simulateMT_default)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset_default)
addParameter(Parser,'centerData',centerData_default)

parse(Parser,varargin{:})

objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
simulateMT = Parser.Results.simulateMT;
theoretical = Parser.Results.theoretical;
equalizeInputsPriorToStimulusOnset = Parser.Results.equalizeInputsPriorToStimulusOnset;
centerData = Parser.Results.centerData;

%%

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
if simulateMT.On
    if isfield(simulateMT,'t')
        t = simulateMT.t;
    else
        t = -500:900;
    end
    speeds = speedsFEF;
    [MT, spref, cohs, MTnull] = simulateMTdata(mtResults,t,simulateMT.modelN,...
        'speedsMT',speeds,'removeBaseline',simulateMT.removeBaseline,'gaussianApprox',simulateMT.gaussianApprox);
    swidth = 1.2*ones(size(spref));
else
    [MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
        'speedsMT',speeds,'cohsMT',cohs,'directionsMT',directionsMT,...
        'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
        'speedPrefOpts',speedPrefOpts);
    
    MTnull = nan(size(MT));
    
    t = mtResults.mt{1}.neuron_t;
end


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);
[interpolatedRnull, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MTnull,speeds,cohs,speedsFEF,cohsFEF);

inputs = permute(interpolatedR,[4,1,2,3]);
inputsNull = permute(interpolatedRnull,[4,1,2,3]);
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
inputsNull = inputsNull(prod(any(isnan(inputsNull),2),[3,4])==0,:,:,:);
sp = spref(~isnan(interpolatedR(1,1,1,:)));%spref(prod(any(isnan(inputs),2),[3,4])==0);
sw = swidth(~isnan(interpolatedR(1,1,1,:)));%swidth(prod(any(isnan(inputs),2),[3,4])==0);

%% Prepare data 
inputs_n = inputs./max(abs(inputs),[],[2,3,4]);
inputsNull_n = inputsNull./max(abs(inputs),[],[2,3,4]);

%% Set MT weights
switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(sp)' - log2(mean(speedsFEF))) ones(size(sp'))];
        
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
            theoreticalInputNull(ri,:,si,ci) = inputsNull_n(:,:,si,ci)'*Atheory(:,ri);
        end
    end
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
theoreticalInputNull = theoreticalInputNull./max(abs(theoreticalInput),[],[2,3,4]);
theoreticalInput = theoreticalInput./max(abs(theoreticalInput),[],[2,3,4]);