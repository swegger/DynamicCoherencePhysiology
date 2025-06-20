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
addParameter(Parser,'cutWindow',[1 1501]);
addParameter(Parser,'t',-100:1600)
addParameter(Parser,'saveType','default')
addParameter(Parser,'acceptInitCohPert',false)

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
acceptInitCohPert = Parser.Results.acceptInitCohPert;

%% Make initCohObject for this file
initCoh = initiateCohObj(dcp.sname,dcp.datapath);
initCoh = assertSpikesExtracted(initCoh,dcp.spikesExtracted);
initCoh = unitsIndex(initCoh);
initCoh = initiateCohTrials(initCoh,trialList,acceptInitCohPert);

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
% initCoh.preferredDirection = directions(maxInds);

if ~isempty(initCoh.unitIndex)
    initCoh = evaluatePreferredDirection(initCoh,win,'speeds',speeds);
    
    % Find rates
    if isnan(speeds)
        spds = unique(initCoh.speeds);
    else
        spds = speeds;
    end
    initCoh.r = calcRates(initCoh,width,'t',t);
    initCoh = cohConditionedRates(initCoh,'width',width,...
        'dirs',NaN,'speeds',spds,...
        'marginalizeDirection',false,'t',t);
    
    
    initCoh.location = dcp.location;
    
    initCoh = findActive(initCoh,rateCutoff,cutWindow);
end

%% Save object to file
if strcmp(initCoh.sname,'ar')
    initCoh.saveLocation = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/'...
        initCoh.datapath(end-8:end-1) 'obj/initCoh' initCoh.datapath(end-8:end)];
elseif strcmp(initCoh.sname,'fr')
    initCoh.saveLocation = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick/'...
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
    