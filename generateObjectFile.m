function dcp = generateObjectFile(sname,baseDir,saveName,varargin)
%% appendObjectFile
%
%   dcp = generateObjectFile(file)
%       Determine data files, run preliminary analysis, makes dcp object
%       list, dcp, and saves the results.
%
%%

%% Defaults
excludeList_default = {'20191021a','20191022a','20191022b','20191111a','20191115a','20191118t','20191119a'};

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'sname')
addRequired(Parser,'baseDir')
addRequired(Parser,'saveName')
addParameter(Parser,'saveFlg',true)
addParameter(Parser,'saveDir',[])
addParameter(Parser,'extractSpikes',true)
addParameter(Parser,'excludeList',excludeList_default)

parse(Parser,sname,baseDir,saveName,varargin{:})

sname = Parser.Results.sname;
baseDir = Parser.Results.baseDir;
saveName = Parser.Results.saveName;
saveFlg = Parser.Results.saveFlg;
saveDir = Parser.Results.saveDir;
extractSpikes = Parser.Results.extractSpikes;
excludeList = Parser.Results.excludeList;

%% Add files
disp('Looking for files to add...')

% Find potential new files
dirOut = dir(baseDir);
regOut = regexpi({dirOut.name},'[0-9]{8}[a-z]{1,3}','match');
ind = 0;
for listi = 1:length(regOut)
    if ~isempty(regOut{listi}) && ~strcmp(regOut{listi}{1}(end-2:end),'plx') && ~strcmp(regOut{listi}{1}(end-1:end),'kk')
        ind = ind+1;
        FileList{ind} = regOut{listi}{1};
    end
end

% Extract data from new files
if exist('FileList','var')
    disp(['Found ' num2str(ind) ' files. Extracting...'])
    dcp = dcpPrelim(sname,FileList,extractSpikes);
else
    disp('No new files found.')
    saveFlg = false;
end


%% Saving
if saveFlg
    if isempty(saveName)
        saveName = datafile;
    end
    
    if isempty(saveDir)
        saveDir = ['~/Projects/DynamicCoherencePhysiology/' sname '/dcpObjects/'];
    end
    
    save([saveDir saveName],'dcp')
else
    disp('dcp object not saved')
end