%% dcpObj
%
%   Defines properties and methods of object used for analysis of
%   DynamicCoherencePhysiology data.
%
%%

classdef dcpObj
    % DynamicCoherencePhysiology analysis class
    properties
        sname;
        datapath;
        spikesExtracted = false;
        trials = 0;
        calib
        unitIndex
        klustaID
        
        trialNumbers;
        trialDataFiles;
        trialtype;
        directions;
        speeds;
        locations;
        eye;
        spikeTimes;
        location;
        
%         dirPref;
%         
%         speedPref;
%         
%         initiateCoh;
%         
%         washout;
%         
%         dynamicCoh;
    end
    
    methods
        
        %% Core methods
        function obj = dcpObj(sname,datapath)
        %Constructor
            obj.sname = sname;
            obj.datapath = datapath;
            
            % For affine transformation of eye signal data
            obj.calib.t = 1:150;
            obj.calib.posGain = 0.025;
            obj.calib.speedGain = 0.09189;
            obj.calib.accThres = 1.1;
            obj.calib.windowSize = 40;
            
        end
        
        function out = returnDatapath(obj)
            out = obj.datapath;
        end
        
        function obj = assertSpikesExtracted(obj,assertion)
            obj.spikesExtracted = assertion;
        end
        
        function obj = extractSpikingData(obj,plxname,plxdir,maestrodir,addSpike)
        % Add spiking data to Maestro data files
            if ~exist('addSpike','var')
                addSpike = false;
            end
            if obj.spikesExtracted
                disp('Spikes already extracted and inserted into Maestro files.')
            else
                if strcmp(plxname(end-2:end),'plx')
                    extractPlxSpikes(plxname,plxdir,maestrodir);
                elseif strcmp(plxname(end-2:end),'pl2')
                    extractPL2Spikes(plxname,plxdir,maestrodir,'addSpike',addSpike);
                else
                    error(['Extraction format ' plxname(end-2:end) ' not recognized!'])
                end
                obj.spikesExtracted = true;
            end
        end
        
        function obj = extractKKData(obj,kkname,kkdir,plxname,plxdir,maestrodir,addSpike)
        % Add spiking data to Maestro data files
            if ~exist('addSpike','var')
                addSpike = false;
            end
            if obj.spikesExtracted
                disp('Spikes already extracted and inserted into Maestro files.')
            else
                if strcmp(plxname(end-2:end),'pl2')
                    [~, ~, ~, ~, unitsIndex] = ...
                        extractKKSpikes(kkname,kkdir,plxname,plxdir,maestrodir,'addSpike',addSpike);
                    obj.klustaID = unitsIndex;
                else
                    error(['Extraction format ' plxname(end-2:end) ' not recognized!'])
                end
                obj.spikesExtracted = true;
            end
        end
        
        function obj = unitsIndex(obj,altName)
            % Allow for different base name
            if ~exist('altName','var')
                altName = obj.sname;
            end
            
            % Find indices in Maestro files that have spike times
            if obj.spikesExtracted
                files = dir(obj.datapath);
                
                % Determine the index of the first data file
                for fileInx = 1:length(files)
                    if length(files(fileInx).name) >= length(altName) && ...
                            strcmp(files(fileInx).name(1:length(altName)),altName)
                        break
                    end
                end
                
                indsList = [];
                for ti = fileInx:length(files)
                    file = readcxdata([obj.datapath '/' files(ti).name]);
                    if iscell(file.sortedSpikes)
                        indsList = [indsList file.sortedSpikes{2}];
%                         indsList = [indsList find(~cellfun(@isempty,file.sortedSpikes))];
                    end
                end
                obj.unitIndex = unique(indsList);
            else
                obj.unitIndex = 1;
            end
        end
        
        function out = myBoxCar(obj,diffs,width)
            out = nan(size(diffs));
            out(abs(diffs) <= width) = 1;
            out(abs(diffs) > width) = 0;
        end
        
        function r = calcRates(obj,width,varargin)
        % Calculate trial by trial rates using box car
            
            % Parse inputs
            Parser = inputParser;
            
            addRequired(Parser,'obj')
            addRequired(Parser,'width')
            addParameter(Parser,'units',NaN)
            addParameter(Parser,'t',-100:1600)
            addParameter(Parser,'trialN',NaN)
            
            parse(Parser,obj,width,varargin{:})
            
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            units = Parser.Results.units;
            t = Parser.Results.t;
            trialN = Parser.Results.trialN;            
            
            if any(isnan(trialN))
                trialN = 1:numel(obj.spikeTimes);
            end
            if any(isnan(units))
                if ~isempty(obj.klustaID)
                    units = obj.klustaID;
                else
                    units = [];
                    for triali = trialN
                        units = [units obj.spikeTimes{triali}{2}];
                    end
                    units = unique(units);
%                     units = 1:numel(obj.unitIndex);
                end
            end
            t = t(:); 
            
            % Find smoothed rates
            f = @(diffs)obj.myBoxCar(diffs,width);
            r = nan(length(t),numel(trialN),numel(units));
            for triali = trialN
                spikeTimes = obj.spikeTimes{triali}(1:2);
                spikeTimes{1} = spikeTimes{1}(ismember(spikeTimes{2},units));
                spikeTimes{2} = spikeTimes{2}(ismember(spikeTimes{2},units));
                [~, ~, ~, rTemp] = spikeTimes2Rate(spikeTimes,...
                    'time',t,'resolution',1,'Filter',f,...
                    'ComputeVariance','Yes','mixedUnits',true,'units',units);
                r(:,triali,:) = rTemp/width;
%                 [~, ~, ~, rTemp] = spikeTimes2Rate(obj.spikeTimes{triali}(units),...
%                     'time',t,'resolution',1,'Filter',f,...
%                     'ComputeVariance','Yes');
%                 r(:,triali,:) = rTemp/width;
            end
        end
        
        function [r,rste] = conditionalRates(obj,width,directions,speeds,t)
            if ~exist('t','var')
                t = -100:1600;
            end
            rAll = obj.calcRates(width,'t',t);
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                sequences,NaN);
            r = rAll(:,condLocial);
            rste = sqrt(var(r,[],2)/sum(condLogical));
        end
        
        function obj = tableImport(obj)
            fileName = [obj.sname obj.datapath(end-6:end)];
            tableFile = [obj.datapath(1:end-9) ...
                'tables/' obj.sname obj.datapath(end-6:end-1) '.xlsx'];
            T = readtable(tableFile);
            indx = find(strcmp(fileName,T.Date));
            obj.location.x = str2num(T.Location_x_y_z_{1}(1:4));
            obj.location.y = str2num(T.Location_x_y_z_{1}(6:9));
            obj.location.z = str2num(T.Location_x_y_z_{1}(11:14));
            obj.location.depth = str2num(T.Location_x_y_z_{indx-1})-...
                str2num(T.Location_x_y_z_{find(strcmp(T.Location_x_y_z_,'Depth'))+1});
        end
        
        %% Analysis methods
        function [condInds, condLogical] = trialSort(obj,directions,speeds,locations,...
                cohs,seqs,perts)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if isnan(directions)
                dMask = true(size(obj.directions));
            else
                dMask = ismember(obj.directions,directions);
            end
            
            if isnan(speeds)
                sMask = true(size(obj.speeds));
            else
                sMask = ismember(obj.speeds,speeds);
            end
            
            if ~exist('locations','var') || any(isnan(locations))
                lMask = true(size(obj.locations));
            else
                for li = 1:size(locations,2)
                    lMask(:,li) = ismember(...
                        obj.locations(:,li),locations(:,li));
                end
            end
            
            if ~exist('cohs','var') || isnan(cohs) 
                cMask = true(size(obj.directions));
            else
                cMask = ismember(obj.coh,cohs);
            end
            
            if ~exist('seqs','var') || isnan(seqs)
                qMask = true(size(obj.directions));
            else
                qMask = ismember(obj.sequences,seqs);
            end
            
            if ~exist('perts','var') || any(isnan(perts))
                pMask = true(size(obj.directions));
            else
                pMask = ismember(obj.perturbations,perts);
            end
            
            condLogical = prod([dMask,sMask,lMask,cMask,qMask,pMask],2,'native');
            condInds = find(condLogical);
        end
                
        function obj = dcpTrials(obj,trialNs,altName)
        % Add spike times from each data file
            if ~exist('altName','var')
                altName = obj.sname;
            end
        
            datapath = obj.returnDatapath;
            files = dir(datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(altName) && ...
                        strcmp(files(fileInx).name(1:length(altName)),altName)
                    break
                end
            end
            
                ind = 0;
                for ti = trialNs
                    
                    if length(files)-fileInx+1 > ti
                        
                        % Read file
                        file = readcxdata([obj.datapath '/' files(ti+fileInx-1).name]);
                        
                        ind = ind+1;
                        % Add spike times
                        if obj.spikesExtracted
                            obj.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);
                        else
                            obj.spikeTimes{ind}{1} = file.spikes;
                        end
                        
                        
                    end
                end
        end
        
        
        function [mE,steE,E] = MeanEyeSpeed(obj,condInds,varargin)
        % Plots mean eye speed for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'condInds')
            addParameter(Parser,'t',NaN)
            addParameter(Parser,'normalizer',1)
            
            parse(Parser,obj,condInds,varargin{:})
            
            obj = Parser.Results.obj;
            condInds = Parser.Results.condInds;
            t = Parser.Results.t;
            normalizer = Parser.Results.normalizer;
            
            E = (sqrt(vertcat(obj.eye(condInds).hvel).^2 + ...
                vertcat(obj.eye(condInds).vvel).^2 ))/normalizer;
            mE = nanmean(E,1);
            steE = sqrt(nanvar(E,[],1)/length(condInds));
            
        end
        
        %% Plotting methods
        function h = rasterPlot(obj,trials,units)
        % Raster plot
%             h = figure;
            if islogical(trials)
                trials = find(trials);
            end
            lineProps.Color = 'k';
            lineProps.LineStyle = '-';
            for triali = 1:length(trials)
                spikeTimes = obj.spikeTimes{trials(triali)}{1};
                spikeTimes = spikeTimes(ismember(obj.spikeTimes{trials(triali)}{2},units));
                if ~isempty(spikeTimes)
                    plotVertical(spikeTimes,...
                        'MinMax',[triali,triali+1],'lineProperties',lineProps);
                end
                hold on
            end
        end
        
        function h = plotMeanEyeVelocity(obj,condLogical,varargin)
        % Plots mean eye velocity for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'condLogical')
            addParameter(Parser,'h',NaN)
            addParameter(Parser,'sh',NaN)
            addParameter(Parser,'t',NaN)
            addParameter(Parser,'color',NaN)
            
            parse(Parser,obj,condLogical,varargin{:})
            
            obj = Parser.Results.obj;
            condLogical = Parser.Results.condLogical;
            h = Parser.Results.h;
            sh = Parser.Results.sh;
            t = Parser.Results.t;
            color = Parser.Results.color;
            
            if ishandle(h)
                figure(h);
            else
                h = figure;
            end
            if ishandle(sh)
                subplot(sh)
            end
            E(:,:,1) = vertcat(obj.eye(~~condLogical).hvel);
            E(:,:,2) = vertcat(obj.eye(~~condLogical).vvel);
            mE = nanmean(E,1);
            steE = sqrt(nanvar(E,[],1)/sum(condLogical));
            
            if isnan(t)
                t = 0:size(mE,2)-1;
            end
            
            patchProps.FaceAlpha = 0.3;
            if ~any(isnan(color))
                patchProps.FaceColor = color;
            else
                color = [0 0 0];
                patchProps.FaceColor = color;
            end
            Etemp = mE(:,:,1);
            steEtemp = steE(:,:,1);
            myPatch(t(:),Etemp(:),steEtemp(:),'patchProperties',patchProps);
            hold on
            Etemp = mE(:,:,2);
            steEtemp = steE(:,:,2);
            myPatch(t(:),Etemp(:),steEtemp(:),'patchProperties',patchProps);
            plot(t,mE(:,:,1),'Color',color,'LineWidth',2)
            plot(t,mE(:,:,2),'--','Color',color,'LineWidth',2)
        end    
        
        function h = plotMeanEyeSpeed(obj,condLogical,varargin)
        % Plots mean eye speed for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'condLogical')
            addParameter(Parser,'h',NaN)
            addParameter(Parser,'sh',NaN)
            addParameter(Parser,'t',NaN)
            addParameter(Parser,'color',NaN)
            addParameter(Parser,'normalizer',1)
            
            parse(Parser,obj,condLogical,varargin{:})
            
            obj = Parser.Results.obj;
            condLogical = Parser.Results.condLogical;
            h = Parser.Results.h;
            sh = Parser.Results.sh;
            t = Parser.Results.t;
            color = Parser.Results.color;
            normalizer = Parser.Results.normalizer;
            
            if ishandle(h)
                figure(h);
            else
                h = figure;
            end
            if ishandle(sh)
                subplot(sh)
            end
            E = (sqrt(vertcat(obj.eye(~~condLogical).hvel).^2 + ...
                vertcat(obj.eye(~~condLogical).vvel).^2 ))/normalizer;
            mE = nanmean(E,1);
            steE = sqrt(nanvar(E,[],1)/sum(condLogical));
            
            if isnan(t)
                t = 0:length(mE)-1;
            end
            
            patchProps.FaceAlpha = 0.3;
            if ~any(isnan(color))
                patchProps.FaceColor = color;
            else
                color = [0 0 0];
                patchProps.FaceColor = color;
            end
            myPatch(t(:),mE(:),steE(:),'patchProperties',patchProps);
            hold on
            plot(t,mE,'Color',color,'LineWidth',2)
        end
        
        
    end
    
end

