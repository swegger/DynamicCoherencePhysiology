%% dynamicCohObj
%
%   Defines properties and methods of dynamic coh trials object used for 
%   analysis of DynamicCoherencePhysiology data.
%
%   Defined as a subclass of the dcpObj
%
%%

classdef dynamicCohObj < dcpObj
    properties
        sequences;
        perturbations;
    end
    
    methods
        %% Core methods
        function obj = dynamicCohObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
        end 
        
        %% dynamicCohTrials
        function obj = dynamicCohTrials(obj,trialNs)
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
                    
                    if length(files)-fileInx+1 > ti
                        
                        % Read file
                        file = readcxdata([obj.datapath '/' files(ti+fileInx-1).name]);
                        
                        trialname = file.trialname;
                        if strcmp(trialname(1:6),'dynCoh')
                            ind = ind+1;
                            % Update trial
                            obj.trialNumbers(ind,1) = ti;
                            obj.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.trialtype(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'q\d{3}');
                            obj.sequences(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'p\d{3}');
                            obj.perturbations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'s\d{3}');
                            obj.speeds(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'x\d{3}');
                            obj.locations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            [startIndex,endIndex] = regexp(trialname,'y\d{3}');
                            obj.locations(ind,2) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            % Add eye information
                            obj.eye(:,ind).hpos = (file.data(1,:) - ...
                                mean(file.data(1,obj.calib.t)))*obj.calib.posGain;
                            obj.eye(:,ind).vpos = (file.data(2,:) - ...
                                mean(file.data(2,obj.calib.t)))*obj.calib.posGain;
                            
                            obj.eye(:,ind).hvel = (file.data(3,:) - ...
                                mean(file.data(3,obj.calib.t)))*obj.calib.speedGain;
                            obj.eye(:,ind).vvel = (file.data(4,:) - ...
                                mean(file.data(4,obj.calib.t)))*obj.calib.speedGain;
                            
                            sacs = saccadeDetect(file.data(3,:)*obj.calib.speedGain,...
                                file.data(4,:)*obj.calib.speedGain,...
                                'accelerationThreshold',obj.calib.accThres,...
                                'windowSize',40);
                            obj.eye(:,ind).hvel(sacs) = NaN;
                            obj.eye(:,ind).vvel(sacs) = NaN;
                            obj.eye(:,ind).saccades = sacs;
                            
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
        end
        
        %% Neural analysis methods
        function [r,rste] = conditionalRates(obj,width,directions,speeds,...
                sequences,perturbations)
            rAll = obj.calcRates(width);
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                sequences,perturbations);
           r = mean(rAll(:,condLogical),2);
           rste = sqrt(var(r,[],2)/sum(condLogical));
        end
        
        function [R,Rste] = dynamicCohSeqConditionedRates(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'width',50)
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'seqs',NaN)
            addParameter(Parser,'perts',NaN)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            dirs = Parser.Results.dirs;
            speeds = Parser.Results.speeds;
            locs = Parser.Results.locs;
            seqs = Parser.Results.seqs;
            perts = Parser.Results.perts;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end            
            if isnan(speeds)
                speeds = unique(obj.speeds);
            end            
%             if any(isnan(locs))
%                 locs = unique(obj.locations,'rows');
%             end            
            if isnan(seqs)
                seqs = unique(obj.sequences);
            end
            if isnan(perts)
                perts = unique(obj.perturbations);
            end
            
            % Condition on sequences
            for seqi = 1:length(seqs)
                [R(:,seqi),Rste(:,seqi)] = conditionalRates(obj,width,...
                    dirs,speeds,seqs(seqi),perts);
            end
            
        end
        
        %% Plotting methods
                
        function h = dynamicCohMeanEyeSpeedDiff(obj,condLogical,controlE,varargin)
        % Plots mean eye speed for a set of trials specified in condLogical
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'condLogical')
            addRequired(Parser,'controlE')
            addParameter(Parser,'h',NaN)
            addParameter(Parser,'sh',NaN)
            addParameter(Parser,'t',NaN)
            addParameter(Parser,'color',NaN)
            
            parse(Parser,obj,condLogical,controlE,varargin{:})
            
            obj = Parser.Results.obj;
            condLogical = Parser.Results.condLogical;
            controlE = Parser.Results.controlE;
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
            stvE = nanvar(E,[],1)/sum(condLogical);
            
            mEc = nanmean(controlE,1);
            stvEc = nanvar(controlE,[],1)/size(controlE,1);
            
            mEd = mE - mEc;
            steEd = sqrt(stvE + stvEc);
            
            if isnan(t)
                t = 0:length(mEd)-1;
            end
            
            patchProps.FaceAlpha = 0.3;
            if ~any(isnan(color))
                patchProps.FaceColor = color;
            else
                color = [0 0 0];
                patchProps.FaceColor = color;
            end
            myPatch(t(:),mEd(:),steEd(:),'patchProperties',patchProps);
            hold on
            plot(t,mEd,'Color',color,'LineWidth',2)
        end
        
        function [h,sh,colors] = dynamicCohMeanEyeSeq(obj,dirs)
        % Plots mean eye speed for each sequence and target speed
            colors = colormap('lines');
            h = gcf;
            speeds = unique(obj.speeds);
            seqs = unique(obj.sequences);
            for si = 1:length(speeds)
                sh(si) = subplot(2,length(speeds),si);
                for seqi = 1:length(seqs)
                    [~,condLogical] = trialSort(obj,dirs,speeds(si),NaN,NaN,seqs(seqi),NaN);
                    plotMeanEyeSpeed(obj,condLogical,'h',h,'sh',sh(si),'color',colors(seqs(seqi),:));
                end
                plotVertical([150 150+0:300:1500]);
                
                sh(si+3) = subplot(2,length(speeds),si+length(speeds));
                [~,condLogicalC] = trialSort(obj,dirs,speeds(si),NaN,NaN,5,NaN);
                controlE = sqrt(vertcat(obj.eye(~~condLogicalC).hvel).^2 + ...
                    vertcat(obj.eye(~~condLogicalC).vvel).^2);
                for seqi = 1:length(seqs)
                    [~,condLogical] = trialSort(obj,dirs,speeds(si),NaN,NaN,seqs(seqi),NaN);
                    dynamicCohMeanEyeSpeedDiff(obj,condLogical,controlE,...
                        'h',h,'sh',sh(si+3),'color',colors(seqs(seqi),:));
                end
                plotVertical([150 150+0:300:1500]);
            end
        end
        
    end
end