%% dcpNeuralAnalysisSummary
%
%
%
%%

%%

%% ar partial correlation
compToBehavior.applyLogrithmicEstimatorCorrection = false;
initiationOrSustained = 'sustained';
switch initiationOrSustained
    case 'initiation'
        gainFitInd = 1;
        tWin = [150 150];
    case 'sustained'
        gainFitInd = 2;
        tWin = [750 750];
end

load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat', 'init')
speeds = [5 10 20]';
slips = -squeeze(init.eye.mean(init.t == 750,:,:)) + speeds;
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
if compToBehavior.applyLogrithmicEstimatorCorrection
    initGain(:,:,2) = initGain(:,:,2).*speeds*0.4./log2(1.4);
    initGain(:,:,3) = initGain(:,:,3).*speeds*0.4./log2(1+0.4.*speeds./slips);
end
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',gainFitInd,'tWin',tWin,...
    'saveFigures',true,'saveResults',true);


%% fr partial correlation
compToBehavior.applyLogrithmicEstimatorCorrection = false;
initiationOrSustained = 'sustained';
switch initiationOrSustained
    case 'initiation'
        gainFitInd = 1;
        tWin = [150 150];
    case 'sustained'
        gainFitInd = 2;
        tWin = [750 750];
end

load('/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat', 'dcp')
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212', 'init')
speeds = [5 10 20]';
slips = -squeeze(init.eye.mean(init.t == 750,:,:)) + speeds;
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
if compToBehavior.applyLogrithmicEstimatorCorrection
    initGain(:,:,2) = initGain(:,:,2).*speeds*0.4./log2(1.4);
    initGain(:,:,3) = initGain(:,:,3).*speeds*0.4./log2(1+0.4.*speeds./slips);
end
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',gainFitInd,'tWin',tWin,...
    'saveFigures',true,'saveResults',true);


%% ar Neuron typing analysis

% To recompute:
% neuronTypingResults('subjects',{'ar'},'resultsFile',false,'initOnly',false,...
%     'saveFigures',true,'saveResults',true)

% To generate results from a saved data file:
neuronTypingResults('subjects',{'ar'},'resultsFile',true,'initOnly',false,...
    'saveFigures',true,'saveResults',false)


%% fr Neuron typing analysis

% To recompute:
% neuronTypingResults('subjects',{'fr'},'resultsFile',false,'initOnly',false,...
%     'saveFigures',true,'saveResults',true)

% To generate results from a saved data file:
neuronTypingResults('subjects',{'fr'},'resultsFile',true,'initOnly',false,...
    'saveFigures',true,'saveResults',false)

%% Both monkeys neuron typing analysis
neuronTypingAnalysis_streamlined({'ar','fr'},...
    'dcpDynCohFile',{[],'dcpObjectsDynCoh20230331'},...
    'chanMap',{[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24],...
    [23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24]},...
    'saveFigures',true,'saveResults',true)

%% ar Neuron typing, initcoh only

% To recompute:
% neuronTypingResults('subjects',{'ar'},'resultsFile',false,'initOnly',true,...
%     'saveFigures',true,'saveResults',true)

% To generate results from a saved data file:
neuronTypingResults('subjects',{'ar'},'resultsFile',true,'initOnly',true,...
    'saveFigures',true,'saveResults',false)


%% fr Neuron typing, initcoh only

% To recompute:
% neuronTypingResults('subjects',{'fr'},'resultsFile',false,'initOnly',true,...
%     'saveFigures',true,'saveResults',true)

% To generate results from a saved data file:
neuronTypingResults('subjects',{'fr'},'resultsFile',true,'initOnly',true,...
    'saveFigures',true,'saveResults',false)

%% Both monkeys neuron typing analysis, initCoh only
neuronTypingAnalysis_streamlined({'ar','fr'},...
    'dcpDynCohFile',{[],'dcpObjectsDynCoh20230331'},...
    'chanMap',{[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24],...
    [23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24]},...
    'initCohCollate',true,'dynCohCollate',false,...
    'saveFigures',true,'saveResults',true)


%% ar targeted dimensionality reduction
plotOpts.On = true;
initCohPertFile = 'initCohPert20240212.mat';          % Upatded analysis removing poor pursuit trials  
neuronTypingComparison.On = true;
neuronTypingComparison.file = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240522.mat';
dimNames = {'Speed','Coherence','Gain','Offset'};

dcpTargetedDimensionalityReduction('dcpObjectFile','dcpObjects20210406',...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
    'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24],...
    'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/' initCohPertFile],...
    'dcpDynCohFile',[],...
    'applyLogrithmicEstimatorCorrection',false,...
    'neuronTypingComparison',neuronTypingComparison,...
    'dimNames',dimNames,...
    'saveFigures',true,...
    'saveResults',true,...
    'plotOpts',plotOpts);

%% fr targeted dimensionality reduction
plotOpts.On = true;
initCohPertFile = 'initCohPert20240212.mat';      % Upatded analysis removing poor pursuit trials
neuronTypingComparison.On = true;
neuronTypingComparison.file = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240522.mat';
dimNames = {'Speed','Coherence','Gain','Offset'};

dcpTargetedDimensionalityReduction('dcpObjectFile','dcpObjects20230322',...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
    'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/' initCohPertFile],...
    'dcpDynCohFile','dcpObjectsDynCoh20230331',...
    'applyLogrithmicEstimatorCorrection',false,...
    'neuronTypingComparison',neuronTypingComparison,...
    'dimNames',dimNames,...
    'saveFigures',true,...
    'saveResults',true,...
    'plotOpts',plotOpts);

%% ar targeted dimensionality reduction, initCoh only
plotOpts.On = true;
initCohPertFile = 'initCohPert20240212.mat';      % Upatded analysis removing poor pursuit trialsneuronTypingComparison.On = true;
neuronTypingComparison.On = true;
neuronTypingComparison.file = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
dimNames = {'Speed','Coherence','Gain','Offset'};

dcpTargetedDimensionalityReduction_initCoh('dcpObjectFile','dcpObjects20210406',...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
    'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
    'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/' initCohPertFile],...
    'applyLogrithmicEstimatorCorrection',false,...
    'optimize2timePoints',false,...
    'neuronTypingComparison',neuronTypingComparison,...
    'dimNames',dimNames,...
    'saveFigures',true,...
    'saveResults',true,...
    'plotOpts',plotOpts);


%% fr targeted dimensionality reduction, initCoh only
plotOpts.On = true;
initCohPertFile = 'initCohPert20240212.mat';      % Upatded analysis removing poor pursuit trialsneuronTypingComparison.On = true;
neuronTypingComparison.file = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
dimNames = {'Speed','Coherence','Gain','Offset'};

dcpTargetedDimensionalityReduction_initCoh('dcpObjectFile','dcpObjects20230322',...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
    'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/' initCohPertFile],...
    'applyLogrithmicEstimatorCorrection',false,...
    'optimize2timePoints',false,...
    'neuronTypingComparison',neuronTypingComparison,...
    'dimNames',dimNames,...
    'saveFigures',true,...
    'saveResults',true,...
    'plotOpts',plotOpts);

%% ar direction preference analysis
dirPrefAnalysis('dcpObjectFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat',...
    'checkUnitType',true)

%% fr direction preference analysis
dirPrefAnalysis('dcpObjectFile','/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat',...
    'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
    'checkUnitType',true)

%% ar MT to FEF regression
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')
speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

trainCondition = [true, false, true;
                  true, false, true;
                  true, false, true];

compToBehavioralGain.On = true;
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat'; %'/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping_initCoh20240212.mat';
compToBehavioralGain.applyLogrithmicEstimatorCorrection = false;

neuronTypingComparison.On = true;
neuronTypingComparison.file = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';

theoretical.weightTheory = 'simple';
theoretical.expansionDef = 'regression';

plotOpts.On = true;

% objectFile = 'allMT_20230419.mat';
objectFile = '~/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';
simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;

MTtoFEFregression(dcp,'speedPrefOpts',speedPrefOpts,'checkMTFit',false,'trainCondition',trainCondition,'rankN',80,...
    'compToBehavioralGain',compToBehavioralGain,...
    'neuronTypingComparison',neuronTypingComparison,'tWin',[-100 100],...
    'objectFile',objectFile,'opponentMT',false,'directionsMT',[0],...
    'theoretical',theoretical,'simulateMT',simulateMT,...
    'plotOpts',plotOpts,'saveResults',true,'saveFigures',true);

%% fr MT to FEF regression
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat', 'dcp')
speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;
trainCondition = [true, false, true;
                  true, false, true;
                  true, false, true];
compToBehavioralGain.On = true;
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212'; %'/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping_initCoh20240212.mat';
compToBehavioralGain.applyLogrithmicEstimatorCorrection = false;

neuronTypingComparison.On = true;
neuronTypingComparison.file = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';

theoretical.weightTheory = 'simple';
theoretical.expansionDef = 'regression';

plotOpts.On = true;

% objectFile = 'allMT_20230419.mat';
objectFile = '~/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';
simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;

MTtoFEFregression(dcp,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'speedPrefOpts',speedPrefOpts,'checkMTFit',false,'trainCondition',trainCondition,'rankN',80,...
    'compToBehavioralGain',compToBehavioralGain,...
    'neuronTypingComparison',neuronTypingComparison,'tWin',[-100 100],...
    'objectFile',objectFile,'opponentMT',false,'directionsMT',[0],...
    'theoretical',theoretical,'simulateMT',simulateMT,...
    'plotOpts',plotOpts,'saveResults',false,'saveFigures',false);

%% ar MT to FEF integrator
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')

speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

trainCondition = [true, true, true;
                  true, true, true;
                  true, true, true];
              
equalizeInputsPriorToStimulusOnset.On = true;
equalizeInputsPriorToStimulusOnset.method = 'time&conditions';
equalizeInputsPriorToStimulusOnset.latency = 30;

compToBehavioralGain.On = true;
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240212.mat';

theoretical.weightTheory = 'simple';
theoretical.expansionDef = 'regression';

mtSampling.resampleN = 10;
mtSampling.sampleN = 81;

modelFEF = fitSimpleFEFmodelToNeurons(dcp,'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 300],'lambdaRidge',10,...
    'opponentMT',false,'directionsMT',0,...
    'tauLB',40,'tauUB',40,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'compToBehavioralGain',compToBehavioralGain,...
    'theoretical',theoretical,...
    'mtSampling',mtSampling);

%% fr MT to FEF integrator
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat')

speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

trainCondition = [true, true, true;
                  true, true, true;
                  true, true, true];
              
equalizeInputsPriorToStimulusOnset.On = true;
equalizeInputsPriorToStimulusOnset.method = 'time&conditions';
equalizeInputsPriorToStimulusOnset.latency = 30;

compToBehavioralGain.On = true;
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20240212.mat';

theoretical.weightTheory = 'optimal';
theoretical.expansionDef = 'bestfit';

mtSampling.resampleN = 10;
mtSampling.sampleN = 81;

modelFEF = fitSimpleFEFmodelToNeurons(dcp,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 120],'lambdaRidge',10,...
    'opponentMT',false,'directionsMT',0,...
    'tauLB',40,'tauUB',40,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'compToBehavioralGain',compToBehavioralGain,...
    'theoretical',theoretical,...
    'mtSampling',mtSampling);

%% ar MT to FEF theoretical model
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')

speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

trainCondition = [true, false, true;
                  true, false, true;
                  true, false, true];
              
equalizeInputsPriorToStimulusOnset.On = true;
equalizeInputsPriorToStimulusOnset.method = 'time&conditions';
equalizeInputsPriorToStimulusOnset.latency = 30;

compToBehavioralGain.On = true;
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240212.mat';

theoretical.weightTheory = 'optimal';
theoretical.expansionDef = 'bestfit';

modelFEF = fitTheoreticalInputModel(dcp,'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 400],...
    'opponentMT',false,'directionsMT',0,...
    'tauLB',40,'tauUB',40,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'compToBehavioralGain',compToBehavioralGain,...
    'theoretical',theoretical,...
    'f',@(x)(x));

%% ar MT to FEF reduced rank model
% load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')

speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

% objectFile = 'allMT_20230419.mat';
objectFile = '~/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';
   
simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;

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

% neuralDataFile =
% '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240212.mat'; % Original file
% neuralDataFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/targetedDimensionality/targetedDimsDynCohAndInitCoh20240529.mat'; % Speed, gain, and offset dimensions only
% neuralDataFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/targetedDimensionality/targetedDimsInitCoh20240610.mat'; % Speed, gain, and offset dimensions only, initCoh only
% temp = load(neuralDataFile,'pInit','initCoh','initGain','dcp','dimNames','meanEyeSpeed','eye_t');
% subject = temp.dcp{1}.sname;
% R = temp.pInit;
% fef_t = temp.initCoh.neuron_t;
% initGain = temp.initGain;
% if isfield(temp,'meanEyeSpeed')
%     meanEyeSpeed = temp.meanEyeSpeed;
%     eye_t = temp.eye_t;
% else
%     meanEyeSpeed = nan(size(R(:,:,:,1)));
%     eye_t = nan(size(R(:,1,1,1)));
% end
% if isfield(temp,'dimNames')
%     dimNames = temp.dimNames;
% else
%     dimNames = {'Speed','Coherence','Gain','Offset'};
% end

neuralDataFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
temp = load(neuralDataFile,'Rinit','initCoh','initGain','subjects','idx','NumClusters','meanEyeSpeed','eye_t','gainRegression');
subject = temp.subjects{1};
for clusi = 1:temp.NumClusters
    Rtemp(:,:,:,clusi) = nanmean(temp.Rinit(:,:,:,temp.idx==clusi),4);
end
R(:,:,:,1) = Rtemp(:,:,:,9);
R(:,:,:,2) = sum(Rtemp.*permute(temp.gainRegression(1).B(1:temp.NumClusters),[4,2,3,1]),4);
fef_t = temp.initCoh.neuron_t;
temp = load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/targetedDimensionality/targetedDimsInitCoh20240610.mat',...
    'initGain','meanEyeSpeed','eye_t');
initGain = temp.initGain;
if isfield(temp,'meanEyeSpeed')
    meanEyeSpeed = temp.meanEyeSpeed;
    eye_t = temp.eye_t;
else
    meanEyeSpeed = nan(size(R(:,:,:,1)));
    eye_t = nan(size(R(:,1,1,1)));
end
dimNames = {'Speed','Gain'};

modelFEF = fitReducedRankModel(subject,R,fef_t,eye_t,initGain,meanEyeSpeed,dimNames,...
    'fitType','fitSigmasExtraInput',...
    'objectFile',objectFile,...
    'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 400],...
    'opponentMT',false,'directionsMT',0,...
    'simulateMT',simulateMT,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'theoretical',theoretical,...
    'centerData',centerData,...
    'saveResults',true,'saveFigures',true);

%% fr MT to FEF reduced rank model
% load('/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat')

speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [128 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [10, 30, 70, 100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

% objectFile = 'allMT_20230419.mat';
objectFile = '~/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';
   
simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;

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

% Load data
neuralDataFile = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20240212.mat';
temp = load(neuralDataFile,'pInit','initCoh','initGain','dcp','dimNames','meanEyeSpeed','eye_t');
subject = temp.dcp{1}.sname;
R = temp.pInit;
fef_t = temp.initCoh.neuron_t;
initGain = temp.initGain;
if isfield(temp,'meanEyeSpeed')
    meanEyeSpeed = temp.meanEyeSpeed;
    eye_t = temp.eye_t;
else
    meanEyeSpeed = nan(size(R(:,:,:,1)));
    eye_t = nan(size(R(:,1,1,1)));
end
if isfield(temp,'dimNames')
    dimNames = temp.dimNames;
else
    dimNames = {'Speed','Coherence','Gain','Offset'};
end

neuralDataFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
temp = load(neuralDataFile,'Rinit','initCoh','initGain','subjects','idx','NumClusters','meanEyeSpeed','eye_t','gainRegression');
subject = temp.subjects{2};
for clusi = 1:temp.NumClusters
    Rtemp(:,:,:,clusi) = nanmean(temp.Rinit(:,:,:,temp.idx==clusi),4);
end
R(:,:,:,1) = Rtemp(:,:,:,9);
R(:,:,:,2) = sum(Rtemp.*permute(temp.gainRegression(2).B(1:temp.NumClusters),[4,2,3,1]),4);
fef_t = temp.initCoh.neuron_t;
temp = load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/targetedDimensionality/targetedDimsInitCoh20240610.mat',...
    'initGain','meanEyeSpeed','eye_t');
initGain = temp.initGain;
if isfield(temp,'meanEyeSpeed')
    meanEyeSpeed = temp.meanEyeSpeed;
    eye_t = temp.eye_t;
else
    meanEyeSpeed = nan(size(R(:,:,:,1)));
    eye_t = nan(size(R(:,1,1,1)));
end
dimNames = {'Speed','Gain'};


modelFEF = fitReducedRankModel(subject,R,fef_t,eye_t,initGain,meanEyeSpeed,dimNames,...
    'fitType','fitSigmas',...
    'objectFile',objectFile,...
    'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 400],...
    'opponentMT',false,'directionsMT',0,...
    'simulateMT',simulateMT,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'theoretical',theoretical,...
    'centerData',centerData,...
    'saveResults',true,'saveFigures',true);

%% eggerDecoder
speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [32 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

compToBehavior.On = true;
compToBehavior.file{1} = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240212.mat';
compToBehavior.file{2} = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20240212.mat';
compToBehavior.applyLogrithmicEstimatorCorrection = false;

theoretical.weightTheory = 'simple';
theoretical.expansionDef = 'bestfit';

simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;

% objectFile = 'allMT_20230419.mat';
objectFile = '~/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';

eggerMTdecoder('objectFile',objectFile,'speedPrefOpts',speedPrefOpts,...
    'compToBehavior',compToBehavior,...
    'theoretical',theoretical,'simulateMT',simulateMT,...
    'estimatorPathwayNonlinearity',@(x)(2.^x),...
    'saveFigures',false,'saveResults',false)