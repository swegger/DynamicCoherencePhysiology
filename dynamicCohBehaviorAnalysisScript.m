%% Defaults

%% Monkey re

% Data from monkey re was taken before formal physiological experiment began. 
% Pass the following variables to analyze data consistent with what was taken from monkey ar 
dcpAccept = true(1,18);
dcpAccept(8:end) = false;
forceRead = true;
pertMap = 1:8;
dynamicCohBehavioralAnalysis('re','dcpObjectsFile','dcpObjects20190923.mat',...
    'forceRead',forceRead,'dcpAccept',dcpAccept,'pertMap',pertMap,'speeds',[8,12,16]);

%% Moneky ar
dyn = dynamicCohBehavioralAnalysis('ar','dcpObjectsFile','dcpObjects20210406.mat');

%% Monkey ar, initiateCoh
[init, gain] = initialCohBehavioralAnalysis('ar','dcpObjectsFile','dcpObjects20210406.mat');