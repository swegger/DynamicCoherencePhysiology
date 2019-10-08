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
        
        trialNumbers;
        trialDataFiles;
        trialtype;
        directions;
        speeds;
        locations;
        eye;
        spikeTimes;
        
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
        
        function obj = extractSpikingData(obj,plxname,plxdir,maestrodir)
        % Add spiking data to Maestro data files
            if obj.spikesExtracted
                disp('Spikes already extracted and inserted into Maestro files.')
            else
                extractPlxSpikes(plxname,plxdir,maestrodir);
                obj.spikesExtracted = true;
            end
        end
        
        function obj = unitsIndex(obj)
        % Find indices in Maestro files that have spike times
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:4),obj.sname)
                    break
                end
            end
            
            indsList = [];
            for ti = fileInx:length(files)
                file = readcxdata([obj.datapath '/' files(fileInx).name]);
                indsList = [indsList find(~cellfun(@isempty,file.sortedSpikes))];
            end
            obj.unitIndex = unique(indsList);
        end
        
        
        %% washoutTrials
        function obj = washoutTrials(obj,trialNs)
        % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:length(obj.sname)),obj.sname)
                    break
                end
            end
            
                ind = 0;
                for ti = trialNs
                    
                    if length(files)-fileInx+1 >= ti
                        
                        % Read file
                        file = readcxdata([obj.datapath '/' files(ti+fileInx-1).name]);
                        
                        trialname = file.trialname;
                        if strcmp(trialname(1:4),'wash')
                            ind = ind+1;
                            % Update trial
                            obj.washout.trialNumbers(ind,1) = ti;
                            obj.washout.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.washout.trialtpe(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'q\d{3}');
                            obj.washout.sequences(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'p\d{3}');
                            obj.washout.perturbations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.washout.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'s\d{3}');
                            obj.washout.speeds(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'x\d{3}');
                            obj.washout.locations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            [startIndex,endIndex] = regexp(trialname,'y\d{3}');
                            obj.washout.locations(ind,2) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            % Add eye information
                            obj.washout.eye(:,ind).hpos = (file.data(1,:) - ...
                                mean(file.data(1,obj.calib.t)))*obj.calib.posGain;
                            obj.washout.eye(:,ind).vpos = (file.data(2,:) - ...
                                mean(file.data(2,obj.calib.t)))*obj.calib.posGain;
                            
                            obj.washout.eye(:,ind).hvel = (file.data(3,:) - ...
                                mean(file.data(3,obj.calib.t)))*obj.calib.speedGain;
                            obj.washout.eye(:,ind).vvel = (file.data(4,:) - ...
                                mean(file.data(4,obj.calib.t)))*obj.calib.speedGain;
                            
                            sacs = saccadeDetect(file.data(3,:)*obj.calib.speedGain,...
                                file.data(4,:)*obj.calib.speedGain,...
                                'accelerationThreshold',obj.calib.accThres,...
                                'windowSize',40);
                            obj.washout.eye(:,ind).hvel(sacs) = NaN;
                            obj.washout.eye(:,ind).vvel(sacs) = NaN;
                            obj.washout.eye(:,ind).saccades = sacs;
                            
                            % Add spike times
                            obj.washout.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex); 
                            if obj.spikesExtracted
                                obj.washout.spikeTimes{ind} = ...
                                    file.sortedSpikes(obj.unitIndex);                                
                            else
                                obj.washout.spikeTimes{ind}{1} = file.spikes;                    
                            end                             
                            
                        end
                        
                    end
                end
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
            
            if ~exist('perts','var') || isnan(perts)
                pMask = true(size(obj.directions));
            else
                pMask = ismember(obj.perturbations,perts);
            end
            
            condLogical = prod([dMask,sMask,lMask,cMask,qMask,pMask],2,'native');
            condInds = find(condLogical);
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
                if ~isempty(obj.spikeTimes{trials(triali)}{units})
                    plotVertical(obj.spikeTimes{trials(triali)}{units},...
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
            E = sqrt(vertcat(obj.eye(~~condLogical).hvel).^2 + ...
                vertcat(obj.eye(~~condLogical).vvel).^2 );
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

