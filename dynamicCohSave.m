function dynamicCohSave(dcp,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'width',50)
addParameter(Parser,'win',[50,350])
addParameter(Parser,'speeds',NaN)
addParameter(Parser,'trialList',1:5000)
addParameter(Parser,'rateCutoff',10);
addParameter(Parser,'saveType','default')

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
width = Parser.Results.width;
win = Parser.Results.win;
speeds = Parser.Results.speeds;
trialList = Parser.Results.trialList;
rateCutoff = Parser.Results.rateCutoff;
saveType = Parser.Results.saveType;

%% Make initCohObject for this file
dynCoh = dynamicCohObj(dcp.sname,dcp.datapath);
dynCoh = assertSpikesExtracted(dynCoh,dcp.spikesExtracted);
dynCoh = unitsIndex(dynCoh);
dynCoh = dynamicCohTrials(dynCoh,trialList);

%% Evaluate preferred direction for each neuron
dirPref = dirPrefObj(dcp.sname,dcp.datapath);
dirPref = assertSpikesExtracted(dirPref,...
    dcp.spikesExtracted);  % Asserts that spiking data has (not) been extracted already
dirPref = unitsIndex(dirPref);                  % Finds the indices to the units
dirPref = dirPrefTrials(dirPref,trialList);        % Finds dirPref trial data
directions = unique(dirPref.directions);
sp = unique(dirPref.speeds);
for di = 1:length(directions)
    counts = conditionalCounts(dirPref,win,directions(di),sp(1));
    countsTotal(:,di) = sum(counts,2);
end
[~,maxInds] = max(countsTotal,[],2);
dynCoh.preferredDirection = directions(maxInds);

directions = unique(dynCoh.directions);
for di = 1:length(directions)
    initCounts = seqConditionedCounts(dynCoh,'dirs',directions(di),'speeds',speeds,'win',win);
    countsTotalInit(:,di) = permute(sum(initCounts,[1,2,3]),[4,1,2,3]);
end
[~,maxInds] = max(countsTotalInit,[],2);
dynCoh.preferredDirectionRelative = directions(maxInds);

%% Find rates
if isnan(speeds)
    spds = unique(dynCoh.speeds);
else
    spds = speeds;
end
dynCoh = dynamicCohSeqConditionedRates(dynCoh,'width',width,...
    'dirs',NaN,'speeds',spds,...
    'marginalizeDirection',false);

dynCoh.r = calcRates(dynCoh,width);

dynCoh.location = dcp.location;

dynCoh.rateCutoff = rateCutoff;

%% Save object to file
if strcmp(dynCoh.sname,'ar')
    dynCoh.saveLocation = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/'...
        dynCoh.datapath(end-8:end-1) 'obj/dynCoh' dynCoh.datapath(end-8:end)];
end
if ~exist(dynCoh.saveLocation(1:end-17),'dir')
    mkdir(dynCoh.saveLocation(1:end-17))
end
switch saveType
    case {'default'}
        save(dynCoh.saveLocation,'dynCoh')
    case {'structure'}
        mysaveobj(obj,dynCoh.saveLocation)     % Saves as stuct; load with myloadobj
end

    