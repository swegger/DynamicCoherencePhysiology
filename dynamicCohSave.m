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
addParameter(Parser,'cutWindow',[1 1501]);
addParameter(Parser,'t',-100:1600)
addParameter(Parser,'saveType','default')

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
width = Parser.Results.width;
win = Parser.Results.win;
speeds = Parser.Results.speeds;
trialList = Parser.Results.trialList;
rateCutoff = Parser.Results.rateCutoff;
cutWindow = Parser.Results.cutWindow;
t = Parser.Results.t;
saveType = Parser.Results.saveType;

%% Make initCohObject for this file
dynCoh = dynamicCohObj(dcp.sname,dcp.datapath);
dynCoh = assertSpikesExtracted(dynCoh,dcp.spikesExtracted);
dynCoh = unitsIndex(dynCoh);
dynCoh = dynamicCohTrials(dynCoh,trialList);
dynCoh = addCoh(dynCoh);

%% Evaluate preferred direction for each neuron
% dirPref = dirPrefObj(dcp.sname,dcp.datapath);
% dirPref = assertSpikesExtracted(dirPref,...
%     dcp.spikesExtracted);  % Asserts that spiking data has (not) been extracted already
% dirPref = unitsIndex(dirPref);                  % Finds the indices to the units
% dirPref = dirPrefTrials(dirPref,trialList);        % Finds dirPref trial data
% directions = unique(dirPref.directions);
% sp = unique(dirPref.speeds);
% for di = 1:length(directions)
%     counts = conditionalCounts(dirPref,win,directions(di),sp(1));
%     countsTotal(:,di) = sum(counts,2);
% end
% [~,maxInds] = max(countsTotal,[],2);
% dynCoh.preferredDirection = directions(maxInds);

if ~isempty(dynCoh.unitIndex)
    dynCoh = evaluatePreferredDirection(dynCoh,win,'speeds',speeds);
    
    % Find rates
    if isnan(speeds)
        spds = unique(dynCoh.speeds);
    else
        spds = speeds;
    end
    dynCoh = dynamicCohSeqConditionedRates(dynCoh,'width',width,...
        'dirs',NaN,'speeds',spds,...
        'marginalizeDirection',false,'t',t);
    
    if isempty(dynCoh.sequences) || length(dynCoh.spikeTimes{1}) < 2
        dynCoh.r = nan(length(dynCoh.neuron_t),length(dynCoh.sequences),length(dynCoh.unitIndex));
    else
        dynCoh.r = calcRates(dynCoh,width,'t',t);
    end
    
    dynCoh.location = dcp.location;
    
    dynCoh = findActive(dynCoh,rateCutoff,cutWindow);
end

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

    