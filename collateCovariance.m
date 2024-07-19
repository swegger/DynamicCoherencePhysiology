function [C,R,N] = collateCovariance(dcp,varargin)
%% collateCovaraince
%
%
%
%%

%% Defaults
interpolation_default.method = 'linear';
interpolation_default.threshold = 100;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'initCohCollate',true)
addParameter(Parser,'dynCohCollate',true)
addParameter(Parser,'winT',[0 1000])
addParameter(Parser,'interpolation',interpolation_default)
addParameter(Parser,'detectPoorPursuitThreshold',1.5)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
initCohCollate = Parser.Results.initCohCollate;
dynCohCollate = Parser.Results.dynCohCollate;
winT = Parser.Results.winT;
interpolation = Parser.Results.interpolation;
detectPoorPursuitThreshold = Parser.Results.detectPoorPursuitThreshold;

%% Collate responses

R = zeros(length(winT(1):winT(2)));
N = 0;

for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
    
        iF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        initCoh = iF.initCoh;
        
        
    if dynCohCollate && initCohCollate        
        % InitCoh data
        iF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        initCoh = iF.initCoh;
        [~,res] = findCovarianceBehavior(initCoh,...
            unique(initCoh.speeds),unique(initCoh.coh),unique(initCoh.directions),...
            winT,interpolation.method,interpolation.threshold,detectPoorPursuitThreshold);
        if ~isempty(res)
            res = res(~any(isnan(res),2),:);
            R = R + res'*res;
            N = N + size(res,1);
        end
        
        % DynCoh data
        dF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        dynCoh = dF.dynCoh;
        [~,res] = findCovarianceBehavior(dynCoh,...
            unique(dynCoh.sequences),unique(dynCoh.perturbations),unique(dynCoh.directions),...
            winT,interpolation.method,interpolation.threshold);
        if ~isempty(res)
            res = res(~any(isnan(res),2),:);
            R = R + res'*res;
            N = N + size(res,1);
        end
        
    elseif dynCohCollate
        % DynCoh data
        dF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        dynCoh = dF.dynCoh;
        [~,res] = findCovarianceBehavior(dynCoh,...
            unique(dynCoh.sequences),unique(dynCoh.perturbations),unique(dynCoh.directions),...
            winT,interpolation.method,interpolation.threshold);
        if ~isempty(res)
            res = res(~any(isnan(res),2),:);
            R = R + res'*res;
            N = N + size(res,1);
        end
        
    elseif initCohCollate
        % InitCoh data
        iF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        initCoh = iF.initCoh;
        [~,res] = findCovarianceBehavior(initCoh,...
            unique(initCoh.speeds),unique(initCoh.coh),unique(initCoh.directions),...
            winT,interpolation.method,interpolation.threshold,detectPoorPursuitThreshold);
        if ~isempty(res)
            res = res(~any(isnan(res),2),:);
            R = R + res'*res;
            N = N + size(res,1);
        end
        
    end
end

%% Calculate covariance
C = R/(N-1);