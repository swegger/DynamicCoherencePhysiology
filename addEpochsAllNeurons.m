function dynamicCoh = addEpochsAllNeurons(dcpfile,dynCohfile,varargin)
%% addEpochsAllNeurons
%
%
%%

%% Defualts
excludeList_default = {'ar/20191021a','ar/20191022a','ar/20191022b'};

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcpfile')
addRequired(Parser,'dynCohfile')
addParameter(Parser,'appendObjects',false)
addParameter(Parser,'trList',1:3500)
addParameter(Parser,'dirs',[0,180])
addParameter(Parser,'saveflg',true)
addParameter(Parser,'excludeList',excludeList_default)
addParameter(Parser,'width',30)

parse(Parser,dcpfile,dynCohfile,varargin{:})

dcpfile = Parser.Results.dcpfile;
dynCohfile = Parser.Results.dynCohfile;
appendObjects = Parser.Results.appendObjects;
trList = Parser.Results.trList;
dirs = Parser.Results.dirs;
saveflg = Parser.Results.saveflg;
excludeList = Parser.Results.excludeList;
width = Parser.Results.width;


%% Add dcp objects
if appendObjects
    % append new dcp objects if unanalyzed files have been found
    dcp = appendObjectFile(dcpfile,'saveName',dcpfile);

else
    % Load dcp object file
    load(dcpfile)
end

%% Compare dcp list and dynCoh list, add epochs to dynCoh objects
if isempty(dynCohfile)
    dynCohfile = ['~/Projects/DynamicCoherencePhysiology/' dcp{1}.sname...
        '/dcpObjects/dynamicCohObjects' datestr(now,'yyyymmdd')];
    for filei = 1:length(dcp)
        dynamicCoh{filei} = dynamicCohObj(dcp{filei}.sname,...
                    dcp{filei}.datapath);
        dynamicCoh{filei} = assertSpikesExtracted(dynamicCoh{filei},...
            dcp{filei}.spikesExtracted);
        dynamicCoh{filei} = unitsIndex(dynamicCoh{filei});
        dynamicCoh{filei} = dynamicCohTrials(dynamicCoh{filei},trList);
        if ~isempty(dynamicCoh{filei}.trialNumbers) && ~any(strcmp(dynamicCoh{filei}.datapath(end-11:end),excludeList))
            disp(['Working on ' dynamicCoh{filei}.datapath(end-11:end)])
            dynamicCoh{filei} = addBehavioralEpochs(dynamicCoh{filei});
            dynamicCoh{filei} = addNeuralEpochs(dynamicCoh{filei},'width',width);
        end
    end
else
    load(dynCohfile)
    for filei = 1:length(dcp)
        for filej = 1:length(dynamicCoh)
            tempflgs(filej) = strcmp(dcp{filei}.datapath,dynamicCoh{filej}.datapath);
        end
        newfile(filei) = ~any(tempflgs);
    end
    newList = find(newfile);
    for filei = 1:length(newList)        
        dynamicCoh{newList(filei)} = dynamicCohObj(dcp{newList(filei)}.sname,...
                    dcp{newList(filei)}.datapath);
        dynamicCoh{newList(filei)} = assertSpikesExtracted(dynamicCoh{newList(filei)},...
            dcp{newList(filei)}.spikesExtracted);
        dynamicCoh{newList(filei)} = unitsIndex(dynamicCoh{newList(filei)});
        dynamicCoh{newList(filei)} = dynamicCohTrials(dynamicCoh{newList(filei)},trList);
        
        if ~isempty(dynamicCoh{newList(filei)}.trialNumbers) && ~any(strcmp(dynamicCoh{newList(filei)}.datapath(end-11:end),excludeList))
            disp(['Working on ' dynamicCoh{newList(filei)}.datapath(end-11:end)])
            dynamicCoh{newList(filei)} = addBehavioralEpochs(dynamicCoh{newList(filei)});
            dynamicCoh{newList(filei)} = addNeuralEpochs(dynamicCoh{newList(filei)},'width',width);
        end
    end
end

%% Save output
if saveflg
    save(dynCohfile,'dynamicCoh')
end