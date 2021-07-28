function initiateCohSave(dcp,varargin)
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
initCoh = initiateCohObj(dcp.sname,dcp.datapath);
initCoh = assertSpikesExtracted(initCoh,dcp.spikesExtracted);
initCoh = unitsIndex(initCoh);
initCoh = initiateCohTrials(initCoh,trialList);

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
initCoh.preferredDirection = directions(maxInds);

directions = unique(initCoh.directions);
for di = 1:length(directions)
    initCounts = cohConditionedCounts(initCoh,'dirs',directions(di),'speeds',speeds,'win',win);
    countsTotalInit(:,di) = permute(sum(initCounts,[1,2,3]),[4,1,2,3]);
end
[~,maxInds] = max(countsTotalInit,[],2);
initCoh.preferredDirectionRelative = directions(maxInds);

%% Find rates
if isnan(speeds)
    spds = unique(initCoh.speeds);
else
    spds = speeds;
end
initCoh.r = calcRates(initCoh,width);
initCoh = cohConditionedRates(initCoh,'width',width,...
    'dirs',NaN,'speeds',spds,...
    'marginalizeDirection',false);


initCoh.location = dcp.location;

initCoh.rateCutoff = rateCutoff;

%% Save object to file
if strcmp(initCoh.sname,'ar')
    initCoh.saveLocation = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/'...
        initCoh.datapath(end-8:end-1) 'obj/initCoh' initCoh.datapath(end-8:end)];
end
if ~exist(initCoh.saveLocation(1:end-17),'dir')
    mkdir(initCoh.saveLocation(1:end-17))
end
switch saveType
    case {'default'}
        save(initCoh.saveLocation,'initCoh')
    case {'structure'}
        mysaveobj(obj,initCoh.saveLocation)     % Saves as stuct; load with myloadobj
end
    