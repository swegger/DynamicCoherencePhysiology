function [Rinit, Rdyn, cellID, passCutoff, locations] = collateFiringRates(dcp,varargin)
%% collateFiringRates
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'dir',[0 180])
addParameter(Parser,'initCohCollate',true)
addParameter(Parser,'dynCohCollate',true)
addParameter(Parser,'preAllocationSize',1000)
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
dir = Parser.Results.dir;
preAllocationSize = Parser.Results.preAllocationSize;
initCohCollate = Parser.Results.initCohCollate;
dynCohCollate = Parser.Results.dynCohCollate;
chanMap = Parser.Results.chanMap;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;

%% Collate responses
passCutoff = nan(preAllocationSize,1);
Rinit = nan(1701,3,3,preAllocationSize);
Rdyn = nan(1701,5,preAllocationSize);
locations = nan(preAllocationSize,3);
cellID = nan(preAllocationSize,100,3);
indx = 1;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
    
    % Add probe info
    dcp{filei} = addProbeInfo(dcp{filei});
    
    if dynCohCollate && initCohCollate
        
        % InitCoh data
        iF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        initCoh = iF.initCoh;
        
        % DynCoh data
        dF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        dynCoh = dF.dynCoh;
        
        if ~isnan(rateCutoff)
            initCoh = findActive(initCoh,rateCutoff,initCoh.cutWindow);
            dynCoh = findActive(dynCoh,rateCutoff,dynCoh.cutWindow);
        end
        if ~isempty(initCoh.R) && (~isempty(dynCoh.R) || ~any(isnan(dynCoh.R(:))) && size(dynCoh.R,2) > 1)
            if length(initCoh.unitIndex) == length(dynCoh.unitIndex)
                
                passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff | dynCoh.passCutoff;
                
                                
                % Get data for each neuron
                if checkUnitType && isprop(dcp{filei},'unitTypes')
                    unitInd = find(strcmp(dcp{filei}.unitTypes,'good'));
                else
                    unitInd = 1:length(initCoh.preferredDirectionRelative);
                end
                for uniti = unitInd
                    ind = find(dir == initCoh.preferredDirectionRelative(uniti));
                    Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
                    
                    ind = find(dir == dynCoh.preferredDirectionRelative(uniti));
                    Rdyn(:,:,indx) = dynCoh.R(:,:,uniti,ind);
                    
                    
                    for j = 1:length(initCoh.unitIndex)
                        cellID(indx,j,1) = filei;
                        cellID(indx,j,2) = initCoh.unitIndex(uniti);
                        cellID(indx,j,3) = initCoh.unitIndex(j);
                    end
                    
                    if isempty(initCoh.location)
                        x = NaN;
                        y = NaN;
                        depth = NaN;
                    elseif length(initCoh.location.x)==24
                        siteIndex = chanMap(initCoh.chansIndex(uniti) == chanMap(:,1),2);
                        x = initCoh.location.x(siteIndex);
                        y = initCoh.location.y(siteIndex);
                        depth = -initCoh.location.depth(siteIndex);
                    elseif length(initCoh.location.x) > 1
                        siteIndex = floor(initCoh.chansIndex(uniti)/4)+1;
                        tempIndex = find(~isnan(initCoh.location.x));
                        if siteIndex>length(tempIndex)
                            x = NaN;
                            y = NaN;
                            depth = NaN;
                        else
                            x = initCoh.location.x(tempIndex(siteIndex));
                            y = initCoh.location.y(tempIndex(siteIndex));
                            depth = -initCoh.location.depth(tempIndex(siteIndex));
                        end
                    else
                        x = initCoh.location.x;
                        y = initCoh.location.y;
                        depth = -initCoh.location.depth;
                    end
                    locations(indx,:) = [x,y,depth];
                    
                    indx = indx+1;
                end
            end
        end
        
        
    elseif dynCohCollate
        % DynCoh data
        dF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        dynCoh = dF.dynCoh;
        
        if ~isnan(rateCutoff)
            dynCoh = findActive(dynCoh,rateCutoff,dynCoh.cutWindow);
        end
        
        if ~isempty(dynCoh.R)
            
            passCutoff(indx:indx+length(dynCoh.passCutoff)-1) = dynCoh.passCutoff;
            
            
            % Get data for each neuron
            if checkUnitType && isprop(dcp{filei},'unitTypes')
                unitInd = find(strcmp(dcp{filei}.unitTypes,'good'));
            else
                unitInd = 1:length(dynCoh.preferredDirectionRelative);
            end
            for uniti = unitInd
                ind = find(dir == dynCoh.preferredDirectionRelative(uniti));
                Rdyn(:,:,indx) = dynCoh.R(:,:,uniti,ind);
                
                for j = 1:length(dynCoh.unitIndex)
                    cellID(indx,j,1) = filei;
                    cellID(indx,j,2) = dynCoh.unitIndex(uniti);
                    cellID(indx,j,3) = dynCoh.unitIndex(j);
                end
                
                if isempty(dynCoh.location)
                    x = NaN;
                    y = NaN;
                    depth = NaN;
                elseif length(dynCoh.location.x)==24
                    siteIndex = chanMap(dynCoh.chansIndex(uniti) == chanMap(:,1),2);
                    x = dynCoh.location.x(siteIndex);
                    y = dynCoh.location.y(siteIndex);
                    depth = -dynCoh.location.depth(siteIndex);
                elseif length(dynCoh.location.x) > 1
                    siteIndex = floor(dynCoh.chansIndex(uniti)/4)+1;
                    tempIndex = find(~isnan(dynCoh.location.x));
                    if siteIndex>length(tempIndex)
                        x = NaN;
                        y = NaN;
                        depth = NaN;
                    else
                        x = dynCoh.location.x(tempIndex(siteIndex));
                        y = dynCoh.location.y(tempIndex(siteIndex));
                        depth = -dynCoh.location.depth(tempIndex(siteIndex));
                    end
                else
                    x = dynCoh.location.x;
                    y = dynCoh.location.y;
                    depth = -dynCoh.location.depth;
                end
                locations(indx,:) = [x,y,depth];
                
                indx = indx+1;
            end
        end
        
        
    elseif initCohCollate
        % InitCoh data
        iF = load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
            dcp{filei}.datapath(end-8:end)]);
        initCoh = iF.initCoh;
        
        if ~isnan(rateCutoff)
            initCoh = findActive(initCoh,rateCutoff,initCoh.cutWindow);
        end
        
        if ~isempty(initCoh.R)
            
            passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;
            
            
            % Get data for each neuron
            if checkUnitType && isprop(dcp{filei},'unitTypes')
                unitInd = find(strcmp(dcp{filei}.unitTypes,'good'));
            else
                unitInd = 1:length(initCoh.preferredDirectionRelative);
            end
            for uniti = unitInd
                ind = find(dir == initCoh.preferredDirectionRelative(uniti));
                Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
                
                for j = 1:length(initCoh.unitIndex)
                    cellID(indx,j,1) = filei;
                    cellID(indx,j,2) = initCoh.unitIndex(uniti);
                    cellID(indx,j,3) = initCoh.unitIndex(j);
                end
                
                if isempty(initCoh.location)
                    x = NaN;
                    y = NaN;
                    depth = NaN;
                elseif length(initCoh.location.x)==24
                    siteIndex = chanMap(initCoh.chansIndex(uniti) == chanMap(:,1),2);
                    x = initCoh.location.x(siteIndex);
                    y = initCoh.location.y(siteIndex);
                    depth = -initCoh.location.depth(siteIndex);
                elseif length(initCoh.location.x) > 1
                    siteIndex = floor(initCoh.chansIndex(uniti)/4)+1;
                    tempIndex = find(~isnan(initCoh.location.x));
                    if siteIndex>length(tempIndex)
                        x = NaN;
                        y = NaN;
                        depth = NaN;
                    else
                        x = initCoh.location.x(tempIndex(siteIndex));
                        y = initCoh.location.y(tempIndex(siteIndex));
                        depth = -initCoh.location.depth(tempIndex(siteIndex));
                    end
                else
                    x = initCoh.location.x;
                    y = initCoh.location.y;
                    depth = -initCoh.location.depth;
                end
                locations(indx,:) = [x,y,depth];
                
                indx = indx+1;
            end
        end    
    end
end

Rinit = Rinit(:,:,:,1:indx-1);
Rdyn = Rdyn(:,:,1:indx-1);
locations = locations(1:indx-1,:);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);
    