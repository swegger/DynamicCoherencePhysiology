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

parse(Parser,varargin{:})

subjects = Parser.Results.subjects;
resultsFile = Parser.Results.resultsFile;
initOnly = Parser.Results.initOnly;

%% Perform typing analysis

for subjecti = 1:length(subjects)
    switch subjects{subjecti}
        case 'ar'
            if resultsFile(subjecti)
                neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20230417.mat')
            else
                %% TODO
            end
        case 'fr'
            if initOnly
                neuronTypingAnalysis_initCoh('dcpObjectFile','dcpObjects20230322',...
                    'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
                    'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
                    'dcpInitCohPertFile','dcpObjectsPert20230817')
            else
                if resultsFile(subjecti)
                    neuronTypingAnalysis('resultsFile','/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohDynCohComp/neuronTyping20230818.mat')
                else
                    neuronTypingAnalysis('dcpObjectFile','dcpObjects20230322',...
                        'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick',...
                        'chanMap',[23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24],...
                        'dcpInitCohPertFile','dcpObjectsPert20230817',...
                        'dcpDynCohFile','dcpObjectsDynCoh20230331')
                end
            end
    end
end