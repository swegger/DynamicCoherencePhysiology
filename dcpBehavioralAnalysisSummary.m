%% dcpBehavioralAnalysisSummary
%
%
%
%%

%% ar initialCoh pertrubation
[init,gain] = initialCohPertBehavioralAnalysis('ar','dcpObjectsFile','dcpObjectsPertTemp',...
    'saveResults',true','saveFigures',false);

%% fr initialCoh perturbation
[init,gain] = initialCohPertBehavioralAnalysis('fr','dcpObjectsFile','dcpObjectsPert20230929',...
    'saveResults',true','saveFigures',false);

%% ar initialCoh

%% fr initialCoh

%% ar dynamicCoh

%% fr dynamicCoh