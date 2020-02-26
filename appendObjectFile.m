function dcp = appendObjectFile(datafile,varargin)
%% appendObjectFile
%
%   dcp = appendObjectFile(file)
%       Determined if there are any unanalyzed data files, runs preliminary
%       analysis, appends results to the dcp object list, dcp, and saves 
%       the results.
%
%%

%% Defaults
excludeList_default = {'20191021a','20191022a','20191022b','20191111a','20191115a','20191118t','20191119a','20200210a','20200213a'};

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'datafile')
addParameter(Parser,'saveName',[])
addParameter(Parser,'saveFlg',true)
addParameter(Parser,'saveDir',[])
addParameter(Parser,'addFile',[])
addParameter(Parser,'extractSpikes',true)
addParameter(Parser,'excludeList',excludeList_default)

parse(Parser,datafile,varargin{:})

datafile = Parser.Results.datafile;
saveName = Parser.Results.saveName;
saveFlg = Parser.Results.saveFlg;
saveDir = Parser.Results.saveDir;
addFile = Parser.Results.addFile;
extractSpikes = Parser.Results.extractSpikes;
excludeList = Parser.Results.excludeList;

%% Load dcp object
load(datafile)

%% Add files
sname = dcp{1}.sname;
baseDir = dcp{1}.datapath(1:end-9);
if isempty(addFile)
    disp('No file specified, looking for files to add...')
    
    % Generate list of files already in dcpObject cell
    for filei = 1:length(dcp)
        nameList{filei} = dcp{filei}.datapath(end-8:end);
    end
        
    % Find potential new files
    dirOut = dir(baseDir);
    regOut = regexpi({dirOut.name},'[0-9]{8}[a-z]{1,3}','match');
    ind = 0;
    for listi = 1:length(regOut)
        if ~isempty(regOut{listi}) && ~strcmp(regOut{listi}{1}(end-2:end),'plx') && ~strcmp(regOut{listi}{1}(end-1:end),'kk')
            ind = ind+1;
            potentialFiles{ind} = regOut{listi}{1};
        end
    end
        
    ind = 0;
    for filei = 1:length(potentialFiles)
        if ~any(strcmp(potentialFiles{filei},nameList)) && ~any(strcmp(potentialFiles{filei},excludeList))
            ind = ind+1;
            FileList{ind} = potentialFiles{filei};
        end
    end
    
    % Extract data from new files
    if exist('FileList','var')
        dcpNew = dcpPrelim(sname,FileList,extractSpikes);
    else
        disp('No new files found.')
    end
else
    % Extract data from desired file
    dcpNew = dcpPrelim(sname,addFile,extractSpikes);
end

% Add to dcpObject cell
if exist('dcpNew','var')
    ind = length(dcp);
    for i = 1:length(dcpNew)
        dcp{ind+i} = dcpNew{i};
    end
else
    saveFlg = false; % If no new files, no reason to save output
end

%% Saving
if saveFlg
    if isempty(saveName)
        saveName = datafile;
    end
    
    if isempty(saveDir)
        saveDir = ['~/Projects/DynamicCoherencePhysiology/' dcp{1}.sname '/dcpObjects/'];
    end
    
    save([saveDir saveName],'dcp')
else
    disp('dcp object not saved')
end