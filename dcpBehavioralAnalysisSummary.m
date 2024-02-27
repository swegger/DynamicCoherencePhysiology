%% dcpBehavioralAnalysisSummary
%
%
%
%%

%% ar initialCoh pertrubation
detectPoorPursuit.On = true;
detectPoorPursuit.threshold = 1.5;
[init,gain] = initialCohPertBehavioralAnalysis('ar','dcpObjectsFile','dcpObjectsPertTemp',...
    'detectPoorPursuit',detectPoorPursuit,'saveResults',true,'saveFigures',true);

%% fr initialCoh perturbation
detectPoorPursuit.On = true;
detectPoorPursuit.threshold = 1.5;
[init,gain] = initialCohPertBehavioralAnalysis('fr','dcpObjectsFile','dcpObjectsPert20230929',...
    'detectPoorPursuit',detectPoorPursuit,'saveResults',true,'saveFigures',true);

%% ar initialCoh

%% fr initialCoh

%% ar dynamicCoh

%% fr dynamicCoh