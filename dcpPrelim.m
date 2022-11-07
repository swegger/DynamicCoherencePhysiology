function dcp = dcpPrelim(subject,FileList,varargin)
%% dcpPrelim
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addRequired(Parser,'FileList')
addParameter(Parser,'extractSpikes',true)
addParameter(Parser,'startChannel','EVT01')
addParameter(Parser,'acceptMU',false)
addParameter(Parser,'checkReceipt',false)

parse(Parser,subject,FileList,varargin{:})

subject = Parser.Results.subject;
FileList = Parser.Results.FileList;
extractSpikes = Parser.Results.extractSpikes;
startChannel = Parser.Results.startChannel;
acceptMU = Parser.Results.acceptMU;
checkReceipt = Parser.Results.checkReceipt;

%% Perform preliminary analysis
for listi = 1:length(FileList)
    data = FileList{listi};
    dataShort = data(3:end);
    datapath = ['/home/seth/Projects/DynamicCoherencePhysiology/' subject '/' data];
    if strcmp(subject,'ar')
        datapath2 = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' data];
    elseif strcmp(subject,'re')
        datapath2 = data;
    elseif strcmp(subject,'fr')
        datapath2 = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick/' data];
    end
    if str2num(dataShort(1:end-1)) > 191114
        plxfile = [subject dataShort '.pl2'];
        kkfile = [subject dataShort '.kwik'];
        kkpath = [datapath2(1:end-1) 'kk'];
        kkflg = exist([kkpath '/' kkfile],'file');
    else
        kkflg = false;
        plxfile = [subject dataShort '.plx'];
    end
    if strcmp(subject,'ar')
        plxpath = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' data(1:end-1) 'plx/'];
    elseif strcmp(subject,'fr')
        plxpath = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick/' data(1:end-1) 'plx/'];
    end
%    plxpath = [datapath(1:end-1) 'plx'];

    dcp{listi} = dcpObj(subject,datapath);
    if extractSpikes && ~strcmp(FileList{listi},'20191021a') && ~strcmp(FileList{listi},'20191022a') && ~strcmp(FileList{listi},'20191022b')
        if kkflg
            dcp{listi} = extractKKData(dcp{listi},kkfile,kkpath,plxfile,plxpath,dcp{listi}.datapath,true,startChannel,acceptMU,checkReceipt);
        else
            dcp{listi} = extractSpikingData(dcp{listi},plxfile,plxpath,dcp{listi}.datapath);
        end
        dcp{listi} = unitsIndex(dcp{listi});                  % Finds the indices to the units
        if kkflg
            dcp{listi} = setUnitType(dcp{listi},kkpath,kkfile);
        end
        dcp{listi} = tableImport(dcp{listi});
    end
end