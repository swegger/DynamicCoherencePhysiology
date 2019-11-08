%% initiateCohObj
%
%   Defines properties and methods of initiate coh trials object used for 
%   analysis of DynamicCoherencePhysiology data.
%
%   Defined as a subclass of the dcpObj
%
%%

classdef initiateCohObj < dcpObj
    properties
        coh;
        cohTuning;
    end
    
    methods
        %% Core methods
        function obj = initiateCohObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
            obj.cohTuning.parameters = nan(1,3);
        end 
        
        %% initiateCohTrials
        function obj = initiateCohTrials(obj,trialNs)
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
                        if strcmp(trialname(1:7),'initCoh')
                            ind = ind+1;
                            % Update triali                            obj.trialNumbers(ind,1) = ti;
                            obj.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.trialtype(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'h\d{3}');
                            obj.coh(ind,1) = ...
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
                cohs)
            rAll = obj.calcRates(width);
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,cohs);
            r = mean(rAll(:,condLogical,:),2);
            rste = sqrt(var(r,[],2)/sum(condLogical));
        end
        
        function [R,Rste] = cohConditionedRates(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'width',50)
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'cohs',NaN)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            dirs = Parser.Results.dirs;
            speeds = Parser.Results.speeds;
            locs = Parser.Results.locs;
            cohs = Parser.Results.cohs;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end            
            if isnan(speeds)
                speeds = unique(obj.speeds);
            end            
%             if any(isnan(locs))
%                 locs = unique(obj.locations,'rows');
%             end            
            if isnan(cohs)
                cohs = unique(obj.coh);
            end
            
            % Condition on sequences
            for speedi = 1:length(speeds)
                for cohi = 1:length(cohs)
                    [R(:,speedi,cohi,:),Rste(:,speedi,cohi,:)] = ...
                        conditionalRates(obj,width,...
                        dirs,speeds(speedi),cohs(cohi));
                end
            end
        end
        
        %% Plotting methods
        function [h,sh,colors] = initiateCohMeanEye(obj,dirs,normalize)
        % Plots mean eye speed for each sequence and target speed
            if ~exist('normalize','var')
                normalize = true;
            end
            h = figure;
            colors = colormap('lines');
            set(h,'Position',[345 557 1965 420]);
            speeds = unique(obj.speeds);
            cohs = unique(obj.coh);
            for hi = 1:length(cohs)
                sh(hi) = subplot(1,length(cohs),hi);
                for si = 1:length(speeds)
                    if normalize
                        nrm = speeds(si);
                    else
                        nrm = 1;
                    end
                    [~,condLogical] = trialSort(obj,dirs,speeds(si),NaN,cohs(hi));
                    plotMeanEyeSpeed(obj,condLogical,'normalizer',nrm,...
                        'h',h,'sh',sh(hi),'color',colors(si,:));
                end
                axis tight
                ax(hi,:) = axis;
            end
            for hi = 1:length(cohs)
                subplot(sh(hi))
                axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
                plotHorizontal(1);
                plotVertical(150);
                xlabel('Time from motion onset (ms)')
                ylabel('$\frac{\textrm{Eye speed}}{\textrm{Target speed}}$')
                mymakeaxis(gca,'xytitle',[num2str(cohs(hi)) '\% \textrm{coherence}'],...
                    'interpreter','latex','yticks',[0.1 1])
            end
            
            for li = 1:(2*length(speeds)+2)
                leg{li} = '';
            end
            for si = 1:length(speeds)
                leg{2*si} = [num2str(speeds(si)) ' deg/s'];
            end
            legend(leg)
        end
    end
end