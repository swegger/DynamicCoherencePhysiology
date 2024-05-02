%% dcpNeuralAnalysisSummary
%
%
%
%%

%%

%% ar partial correlation
compToBehavior.applyLogrithmicEstimatorCorrection = true;

load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat', 'init')
speeds = [5 10 20]';
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
if compToBehavior.applyLogrithmicEstimatorCorrection
    initGain(:,:,2) = initGain(:,:,2).*speeds*0.4./log2(1.4);
    initGain(:,:,3) = initGain(:,:,3).*speeds*0.4./log2(0.4*speeds);
end
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',1,'tWin',[150 150],...
    'saveFigures',true,'saveResults',true);


%% fr partial correlation
compToBehavior.applyLogrithmicEstimatorCorrection = true;

load('/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat', 'dcp')
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212', 'init')
speeds = [5 10 20]';
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
if compToBehavior.applyLogrithmicEstimatorCorrection
    initGain(:,:,2) = initGain(:,:,2).*speeds*0.4./log2(1.4);
    initGain(:,:,3) = initGain(:,:,3).*speeds*0.4./log2(0.4*speeds);
end
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',1,'tWin',[150 150],...
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
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping_initCoh20240212.mat';
compToBehavior.applyLogrithmicEstimatorCorrection = true;

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
    'compToBehavioralGain',compToBehavioralGain,'tWin',[-100 100],...
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
compToBehavioralGain.On = false;
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping_initCoh20240212.mat';
compToBehavior.applyLogrithmicEstimatorCorrection = true;

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
    'compToBehavioralGain',compToBehavioralGain,'tWin',[-100 100],...
    'objectFile',objectFile,'opponentMT',false,'directionsMT',[0],...
    'theoretical',theoretical,'simulateMT',simulateMT,...
    'plotOpts',plotOpts,'saveResults',true,'saveFigures',true);

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
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')

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

modelFEF = fitReducedRankModel('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240212.mat',...
    'objectFile',objectFile,...
    'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 400],...
    'opponentMT',false,'directionsMT',0,...
    'simulateMT',simulateMT,...
    'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
    'theoretical',theoretical,...
    'centerData',centerData);

%% eggerDecoder
speedPrefOpts.tWin = [40,120];
speedPrefOpts.P0 = [16,1];
speedPrefOpts.ub = [32 128];
speedPrefOpts.lb = [1,0];
speedPrefOpts.c = [100];
speedPrefOpts.s = NaN;
speedPrefOpts.d = 0;

compToBehavior.On = true;
compToBehavior.file{1} = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPertMeansGains20240212';
compToBehavior.file{2} = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPertMeansGains20240212';
compToBehavior.applyLogrithmicEstimatorCorrection = true;

theoretical.weightTheory = 'simple';
theoretical.expansionDef = 'bestfit';

simulateMT.On = true;
simulateMT.modelN = 1000;
simulateMT.removeBaseline = true;
simulateMT.gaussianApprox = true;

% objectFile = 'allMT_20230419.mat';
objectFile = '~/Projects/DynamicCoherencePhysiology/ar/mtdata/temp.mat';

eggerMTdecoder('objectFile',objectFile,'speedPrefOpts',speedPrefOpts,'compToBehavior',compToBehavior,...
    'theoretical',theoretical,'simulateMT',simulateMT,'saveFigures',true,'saveResults',true)