function neuronTypingResults(varargin)
%% neuronTypingResults
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'subjects',{'ar','fr'})
addParameter(Parser,'resultsFile',[true false])
addParameter(Parser,'initOnly',false)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,varargin{:})

subjects = Parser.Results.subjects;
resultsFile = Parser.Results.resultsFile;
initOnly = Parser.Results.initOnly;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Perform typing analysis

for subjecti = 1:length(subjects)
    switch subjects{subjecti}
        case 'ar'
            if initOnly
                if resultsFile(subjecti)
                    
                    % Updated analysis
                    neuronTypingAnalysis_initCoh('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping_initCoh20240209.mat',...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);  
                    
                    % Updated analysis, removing poor pursuit trials from
                    neuronTypingAnalysis_initCoh('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping_initCoh20240212.mat',...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);                    
                else
%                     initCohPertFile = initCohPertTemp.mat;            % Original result of initialCohPertBehavioralAnalaysis
%                     initCohPertFile = initCohPert20240208.mat;        % Updated analysis initialCohPertBehavioralAnalaysis
                    initCohPertFile = initCohPert20240212.mat;          % Upatded analysis removing poor pursuit trials  
                    neuronTypingAnalysis_initCoh('dcpObjectFile','dcpObjects20210406',...
                        'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
                        'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24],...
                        'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/' initCohPertFile],...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);
                end
            else
                
                if resultsFile(subjecti)
                    % Original analysis
%                     neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20230417.mat',...
%                         'saveFigures',saveFigures,...
%                         'saveResults',false);

                    % Updated analysis
%                     neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240208.mat',...
%                         'saveFigures',saveFigures,...
%                         'saveResults',saveResults);

                    % Updated analysis, removing poor pursuit trials from
                    neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20240212.mat',...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);
                else   
%                     initCohPertFile = initCohPertTemp.mat;            % Original result of initialCohPertBehavioralAnalaysis
%                     initCohPertFile = initCohPert20240208.mat;        % Updated analysis initialCohPertBehavioralAnalaysis
                    initCohPertFile = 'initCohPert20240212.mat';          % Upatded analysis removing poor pursuit trials                 
                    neuronTypingAnalysis('dcpObjectFile','dcpObjects20210406',...
                        'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle',...
                        'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24],...
                        'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/' initCohPertFile],...
                        'includeEarlyInitCohPertTime',true,...
                        'dcpDynCohFile',[],...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);
                end
            end
        case 'fr'
            if initOnly
                if resultsFile(subjecti)
                    neuronTypingAnalysis_initCoh('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping_initCoh20240212.mat',...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);                    
                else
                    neuronTypingAnalysis_initCoh('dcpObjectFile','dcpObjects20230322',...
                        'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
                        'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
                        'initCohPertFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat',...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);
                end
            else
                if resultsFile(subjecti)
                    % Original analysis
%                     neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20230818.mat',...
%                         'saveFigures',saveFigures,...
%                         'saveResults',false);

                    % Updated analysis
%                     neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20240208.mat',...
%                         'saveFigures',saveFigures,...
%                         'saveResults',saveResults);

                    % Updated analysis, removing poor pursuit trials from
                    % estimate of gain during initiate coh
                    neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20240212.mat',...
                        'includeEarlyInitCohPertTime',true,...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);
                else
%                     initCohPertFile = initCohPert20230822.mat;      % Original result of initialCohPertBehavioralAnalaysis
%                     initCohPertFile = initCohPert20240208.mat;      % Updated analysis initialCohPertBehavioralAnalaysis
                    initCohPertFile = 'initCohPert20240212.mat';      % Upatded analysis removing poor pursuit trials
                    neuronTypingAnalysis('dcpObjectFile','dcpObjects20230322',...
                        'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
                        'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
                        'initCohPertFile',['/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/' initCohPertFile],...
                        'includeEarlyInitCohPertTime',true,...
                        'dcpDynCohFile','dcpObjectsDynCoh20230331',...
                        'saveFigures',saveFigures,...
                        'saveResults',saveResults);
                end
            end
    end
end