function predictDynCohRRdynamicsModel(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'neuronTypingFile',NaN)
addParameter(Parser,'fitReducedRankDynamicsFile',NaN)
addParameter(Parser,'MTtoFEFregressionFile',NaN)
addParameter(Parser,'speedInd',2)
addParameter(Parser,'cohInds',[2 1 3])
addParameter(Parser,'transitionTime',490)
addParameter(Parser,'timeconstant',20)
addParameter(Parser,'inputGain',[10 2])

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
neuronTypingFile = Parser.Results.neuronTypingFile;
fitReducedRankDynamicsFile = Parser.Results.fitReducedRankDynamicsFile;
MTtoFEFregressionFile = Parser.Results.MTtoFEFregressionFile;
speedInd = Parser.Results.speedInd;
cohInds = Parser.Results.cohInds;
transitionTime = Parser.Results.transitionTime;
timeconstant = Parser.Results.timeconstant;
inputGain = Parser.Results.inputGain;


%% Check for specific files as input, otherwise set default files
if any(isnan(neuronTypingFile))
    switch subject
        case 'ar'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
        case 'fr'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';            
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(fitReducedRankDynamicsFile))
    switch subject
        case 'ar'
%             fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240613.mat';
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240729.mat';
        case 'fr'
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/fitReducedRankModel20240613.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(MTtoFEFregressionFile))
    switch subject
        case 'ar'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/ar/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        case 'fr'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/fr/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

%% Load the data
neuronTyping = load(neuronTypingFile,'subjects','Y','cellID','gainRegression','idx','Rinit','NumClusters','initCoh');
rrModel = load(fitReducedRankDynamicsFile,'modelFEF','theoreticalInput','kappas');
MTtoFEF = load(MTtoFEFregressionFile,'vOpt');
t = rrModel.modelFEF.t;
modelFEF = rrModel.modelFEF;

%% Create weighted MT response for 60% coh to transitions to 20 and 100% coh
inputInitCoh = rrModel.theoreticalInput;
input = inputInitCoh(:,:,speedInd,cohInds(1));
transitionInd = find(t == transitionTime);
for i = 2:length(cohInds)
    input(:,:,i) = inputStreamInterpolate(inputInitCoh(:,:,speedInd,cohInds(1)),...
        inputInitCoh(:,:,speedInd,cohInds(i)),...
        transitionInd,timeconstant,inputGain);
end
input(3,:,:) = inputInitCoh(3,:,3,1:3);

%% Run simulation of DynCoh
for i = 1:length(cohInds)
    [~,kappas(:,:,i)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
        't',modelFEF.t,'us',input(:,:,i),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
end

%% Plot
figure
for ki = 1:size(kappas,1)
    subplot(1,size(kappas,1),ki)
    for ci = 1:size(rrModel.kappas,4)
        plot(t,rrModel.kappas(ki,:,speedInd,ci),'k')
        hold on
    end
    for ci = 2:length(cohInds)
        plot(t,kappas(ki,:,ci))
    end
    xlabel('Time from motion onset (ms)')
end

figure
for inputi = 1:size(input,1)
    subplot(1,size(input,1),inputi)
    plot(t,squeeze(inputInitCoh(inputi,:,speedInd,cohInds)),'k')
    hold on
    ind = 1:size(inputInitCoh,3);
    plot(t,squeeze(inputInitCoh(inputi,:,ind~=speedInd,3)),'k--')
    plot(t,squeeze(inputInitCoh(inputi,:,ind~=speedInd,2)),'k--')
    plot(t,squeeze(inputInitCoh(inputi,:,ind~=speedInd,1)),'k--')
    plot(t,squeeze(input(inputi,:,2:end)))
end


%% Functions
function input = inputStreamInterpolate(stream1,stream2,transitionInd,timeconstant,inputGain)
    %%
    t = 0:100;
    f = exp(-t/timeconstant);
    f = f/sum(f);
    input = stream1;
    streamDiff = stream2(:,transitionInd:end)-stream1(:,transitionInd:end);
    for inputi = 1:size(input,1)
        streamDiff_pos = zeros(size(streamDiff(inputi,:)));
        streamDiff_pos(streamDiff(inputi,:) > 0) = streamDiff(inputi,streamDiff(inputi,:) > 0);
        streamDiff_neg = zeros(size(streamDiff(inputi,:)));
        streamDiff_neg(streamDiff(inputi,:) < 0) = streamDiff(inputi,streamDiff(inputi,:) < 0);
        input(inputi,transitionInd:end) = input(inputi,transitionInd:end) + ...
            inputGain(1)*streamDiff_pos + inputGain(2)*streamDiff_neg;
        %     input(:,transitionInd-length(t):end) = filter(f,1,input(:,transitionInd-length(t):end),[],2);
        input(inputi,:) = lowpass(input(inputi,:),1/timeconstant,1000);
    end