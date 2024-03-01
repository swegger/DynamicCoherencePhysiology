%% dcpNeuralAnalysisSummary
%
%
%
%%

%%

%% ar partial correlation
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat', 'init')
speeds = [5 10 20]';
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',2,'saveFigures',true);


%% fr partial correlation
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat', 'dcp')
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212', 'init')
speeds = [5 10 20]';
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',2,'saveFigures',true);


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
theoretical.weightTheory = 'optimal';
theoretical.expansionDef = 'bestfit';
plotOpts.On = true;
MTtoFEFregression(dcp,'speedPrefOpts',speedPrefOpts,'checkMTFit',false,'trainCondition',trainCondition,'rankN',80,...
    'compToBehavioralGain',compToBehavioralGain,'tWIn',[-100 100],...
    'opponentMT',false,'directionsMT',[0],...
    'theoretical',theoretical,...
    'plotOpts',plotOpts);

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
compToBehavioralGain.file = '/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping_initCoh20240212.mat';
theoretical.weightTheory = 'optimal';
theoretical.expansionDef = 'bestfit';
plotOpts.On = true;
MTtoFEFregression(dcp,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
    'speedPrefOpts',speedPrefOpts,'checkMTFit',false,'trainCondition',trainCondition,'rankN',80,...
    'compToBehavioralGain',compToBehavioralGain,'tWIn',[-100 100],...
    'opponentMT',false,'directionsMT',[0],...
    'theoretical',theoretical,...
    'plotOpts',plotOpts);

%% ar MT to FEF integrator
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
modelFEF = fitSimpleFEFmodelToNeurons(dcp,'speedPrefOpts',speedPrefOpts,...
    'checkMTFit',false,'trainCondition',trainCondition,...
    'tWIn',[-100 100],...
    'opponentMT',false,'directionsMT',0);