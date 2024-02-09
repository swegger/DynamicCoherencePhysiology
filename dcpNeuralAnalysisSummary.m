%% dcpNeuralAnalysisSummary
%
%
%
%%

%%

%% ar partial correlation
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat')
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240208.mat', 'init')
speeds = [5 10 20]';
initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
clear init
[RHO, PVAL, Z] = neuronPartialCorrelation(dcp,...
    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
    'regressionModels',{'speed','gain'},'gains',initGain(:,:,2:3),...
    'tGainMeasurement',[150 750],'gainFitInd',2,'saveFigures',true);


%% fr partial correlation
load('dcpObjects20230322.mat', 'dcp')
load('initCohPert20240208', 'init')
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