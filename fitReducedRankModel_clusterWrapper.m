function fitReducedRankModel_clusterWrapper(subject,fitType,fitNumber,varargin)
%%
%
%
%
%
%%

%% Defaults
neuralDataFile_default = '/hpc/group/lisbergerlab/se138/Projects/DynamicCoherencePhysiology/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
plotOpts_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addRequired(Parser,'fitType')
addRequired(Parser,'fitNumber')
addParameter(Parser,'neuralDataFile',neuralDataFile_default)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,subject,fitType,fitNumber,varargin{:})

subject = Parser.Results.subject;
fitType = Parser.Results.fitType;
fitNumber = Parser.Results.fitNumber;
neuralDataFile = Parser.Results.neuralDataFile;
plotOpts = Parser.Results.plotOpts;

%% Set file names and locations
behavioralDataFile= ['/hpc/group/lisbergerlab/se138/Projects/DynamicCoherencePhysiology/'...
    subject '/targetedDimensionality/targetedDimsInitCoh20240610.mat'];

saveFileFull = ['/hpc/group/lisbergerlab/se138/Projects/DynamicCoherencePhysiology/'...
    subject '/ReducedRankModelFits/' fitType '/fitReducedRankModel_' fitType '_' num2str(fitNumber) '_' datestr(now,'yyyymmdd')];
        

%% Set up the variabiles

objectFile = '/hpc/group/lisbergerlab/se138/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';

speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;
simulateMT.accept = true(100,1);
simulateMT.accept(end-26:end) = false;

equalizeInputsPriorToStimulusOnset.On = false;
equalizeInputsPriorToStimulusOnset.method = 'time&conditions';
equalizeInputsPriorToStimulusOnset.latency = 30;

theoretical.weightTheory = 'simple';
theoretical.expansionDef = 'bestfit';

% Train condition not yet implemented!
trainCondition = [true, false, true;
                  true, false, true;
                  true, false, true];
              
centerData.On = false;
centerData.inds = [2,2];

dimNames = {'Speed','Gain'};

% Load data
temp = load(neuralDataFile,'Rinit','initCoh','initGain','idx','NumClusters','meanEyeSpeed','eye_t','gainRegression');
for clusi = 1:temp.NumClusters
    Rtemp(:,:,:,clusi) = nanmean(temp.Rinit(:,:,:,temp.idx==clusi),4);
end
R(:,:,:,1) = Rtemp(:,:,:,9);
R(:,:,:,2) = sum(Rtemp.*permute(temp.gainRegression(1).B(1:temp.NumClusters),[4,2,3,1]),4);
fef_t = temp.initCoh.neuron_t;
temp = load(behavioralDataFile,'initGain','meanEyeSpeed','eye_t');
initGain = temp.initGain;
if isfield(temp,'meanEyeSpeed')
    meanEyeSpeed = temp.meanEyeSpeed;
    eye_t = temp.eye_t;
else
    meanEyeSpeed = nan(size(R(:,:,:,1)));
    eye_t = nan(size(R(:,1,1,1)));
end

loadFromDataFile.On = false;


%% Fit the model
modelFEF = fitReducedRankModel(subject,R,fef_t,eye_t,initGain,meanEyeSpeed,dimNames,...
    'fitType',fitType,...
    'objectFile',objectFile,...
    'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 400],...
    'opponentMT',false,'directionsMT',0,...
    'simulateMT',simulateMT,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'theoretical',theoretical,...
    'centerData',centerData,...
    'saveResults',true,'saveLocation',saveFileFull,...
    'saveFigures',false,'loadFromDataFile',loadFromDataFile,...
    'plotOpts',plotOpts);
