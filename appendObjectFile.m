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
excludeList_default = {'20191021a','20191022a','20191022b','20191111a',...
    '20191115a','20191118t','20191119a','20200210a','20200213a','20200218a',...
    '20200221a','20200309b','20200316a','20210205a','20200214a','20200214b',...
    '20210203a','20210308a'};

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'datafile')
addParameter(Parser,'saveName',[])
addParameter(Parser,'saveFlg',true)
addParameter(Parser,'saveDir',[])
addParameter(Parser,'addFile',[])
addParameter(Parser,'extractSpikes',true)
addParameter(Parser,'excludeList',excludeList_default)
addParameter(Parser,'assertkk',false)
addParameter(Parser,'startChannel','EVT01')
addParameter(Parser,'acceptMU',false)
addParameter(Parser,'checkReceipt',false)

parse(Parser,datafile,varargin{:})

datafile = Parser.Results.datafile;
saveName = Parser.Results.saveName;
saveFlg = Parser.Results.saveFlg;
saveDir = Parser.Results.saveDir;
addFile = Parser.Results.addFile;
extractSpikes = Parser.Results.extractSpikes;
excludeList = Parser.Results.excludeList;
assertkk = Parser.Results.assertkk;
startChannel = Parser.Results.startChannel;
acceptMU = Parser.Results.acceptMU;
checkReceipt = Parser.Results.checkReceipt;

%% Load dcp object
load(datafile)

%% Add files
sname = dcp{1}.sname;
baseDir = [dcp{1}.datapath(1:end-12) sname '/'];
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
    
    % Check for kwik files if assertkk
    if assertkk
        for filei = 1:length(FileList)
            data = FileList{filei};
            dataShort = data(3:end);
            if strcmp(sname,'ar')
                datapath2 = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' data];
            end
            kkfile = [sname dataShort '.kwik'];
            kkpath = [datapath2(1:end-1) 'kk'];
            kkflg(filei) = exist([kkpath '/' kkfile],'file');
        end
        FileList = FileList(kkflg);
    end
    
    % Extract data from new files
    if exist('FileList','var')
        ind = length(dcp);
        for filei = 1:length(FileList)
            disp(['Working on ' FileList{filei} ', file ' num2str(filei) ' of ' num2str(length(FileList)) '.'])
            dcpNew = dcpPrelim(sname,FileList(filei),'extractSpikes',extractSpikes,'startChannel',startChannel,'acceptMU',acceptMU,'checkReceipt',checkReceipt);
            dcp{ind+filei} = dcpNew{1};
            
            % Saving
            if saveFlg
                if isempty(saveName)
                    saveName = datafile;
                end
                
                if isempty(saveDir)
                    saveDir = ['~/Projects/DynamicCoherencePhysiology/' dcp{1}.sname '/dcpObjects/'];
                end
                
                save([saveDir saveName],'dcp')
                disp(['File ' FileList{filei} ' added and saved to ' [saveDir saveName] ' object list.'])
            else
                disp('dcp object not saved')
            end
        end
    else
        disp('No new files found.')
    end
else
    % Extract data from desired file
    dcpNew = dcpPrelim(sname,addFile,'extractSpikes',extractSpikes,'startChannel',startChannel,'acceptMU',acceptMU,'checkReceipt',checkReceipt);
    
    
    % Add to dcpObject cell
    if exist('dcpNew','var')
        ind = length(dcp);
        for i = 1:length(dcpNew)
            dcp{ind+i} = dcpNew{i};
        end
    else
        saveFlg = false; % If no new files, no reason to save output
    end
    
    % Saving
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
end


