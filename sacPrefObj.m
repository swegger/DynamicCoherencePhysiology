%% sacPrefObj
%
%   Defines properties and methods of saccade direction preference object 
%   used for analysis of DynamicCoherencePhysiology data.
%
%   Defined as a subclass of the dcpObj
%
%%

classdef sacPrefObj < dcpObj
    properties (SetAccess = private)
        objType = 'sacPrefObj';
        directionTuning;
        r;
        reward;
        rotation;
        amplitudes;
        
        filterWidth;
        R;
        Rste;
        eye_t;
        neuron_t;
        rateCutoff;
        cutWindow;
        passCutoff;
        
    end
    
    methods
        %% Core methods
        function obj = sacPrefObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
            obj.directionTuning.tuningFun = @vonMises;
            obj.directionTuning.parameters = nan(1,3);
        end  
        
        %% sacPrefTrials
        function obj = sacPrefTrials(obj,trialNs)
        % Add direction preference trials to data object
            datapath = obj.returnDatapath;
            files = dir(datapath);
            
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
                        if ~isempty(trialname) && strcmp(trialname(1:7),'sacPref')
                            ind = ind+1;
                            % Update trial
                            obj.trialNumbers(ind,1) = ti;
                            obj.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.trialtype(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'a\d{3}');
                            obj.amplitudes(ind,1) = ...
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
                            obj.eye(:,ind).saccades = sacs;
                            
                            obj.rotation(ind,1) = file.key.iVelTheta/1000;
                            
                            onTemp = [file.targets.on{2:end}];
                            t_temp = min(onTemp):(max(onTemp)-1);
                            obj.eye(:,ind).t = t_temp;
                            acquiredInd = find(sqrt( ...
                                (obj.eye(:,ind).hpos-obj.amplitudes(ind,1)*cosd(obj.directions(ind,1)+obj.rotation(ind,1))).^2 + ...
                                (obj.eye(:,ind).vpos-obj.amplitudes(ind,1)*sind(obj.directions(ind,1)+obj.rotation(ind,1))).^2) < 2,1);
                            if ~isempty(acquiredInd)
                                obj.eye(:,ind).targetAcquired = t_temp(acquiredInd);
                            else
                                obj.eye(:,ind).targetAcquired = NaN;
                            end
                            
                            % Add other information
                            obj.reward(ind,1) = file.key.iRewLen1;
                            
                            % Add spike times
                            if obj.spikesExtracted
                                obj.spikeTimes{ind} = ...
                                    file.sortedSpikes;  
%                                 obj.spikeTimes{ind} = ...
%                                     file.sortedSpikes(obj.unitIndex);                                
                            else
                                obj.spikeTimes{ind}{1} = file.spikes;     
                                obj.spikeTimes{ind}{2} = ones(size(file.spikes));
                            end                         
                            
                        end
                        
                    end
                end
        end
        
        %% Neural analysis methods
        function [r,rste] = conditionalRates(obj,width,directions,t,alignToSaccade)
            if ~exist('t','var')
                t = -1000:100;
            end
            if ~exist('alignToSaccade')
                alignToSaccade = true;
                t_offsets = [obj.eye(:).targetAcquired];
            end
            if alignToSaccade
                rAll = obj.calcRates(width,'t',t,'t_offsets',t_offsets);
            else
                rAll = obj.calcRates(width,'t',t);
            end
            [~,condLogical] = trialSort(obj,directions,NaN,NaN);
            r = mean(rAll(:,condLogical,:),2);
            rste = sqrt(var(rAll(:,condLogical,:),[],2)/sum(condLogical));
        end
        
        function obj = dirConditionedRates(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'width',50)
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'t',-1000:100)
            addParameter(Parser,'alignToSaccade',true)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            dirs = Parser.Results.dirs;
            t = Parser.Results.t;
            alignToSaccade = Parser.Results.alignToSaccade;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end    
            
            % Condition on sequences
            for di = 1:length(dirs)
                [R(:,di,:),Rste(:,di,:)] = ...
                    conditionalRates(obj,width,...
                    dirs(di),t,alignToSaccade);
            end
            if exist('R','var')
                obj.R = R;
                obj.Rste = Rste;
            else
                obj.R = [];
                obj.Rste = [];
            end
            obj.filterWidth = width;
            obj.neuron_t = t;
            
        end
        
        function obj = dirRates(obj,boxCarWidth,t_offsets)
            if exist('t_offsets','var')
                obj.r = obj.calcRates(boxCarWidth,'t',-1000:100,'t_offsets',t_offsets);
            else
                obj.r = obj.calcRates(boxCarWidth);
            end
        end
        
        
        function counts = conditionalCounts(obj,win,directions)
            countsAll = obj.countSpikes(win);
            if exist('countsAll','var')
                [~,condLogical] = trialSort(obj,directions,NaN,NaN,NaN);
                counts = countsAll(:,condLogical);
            else
                [~,condLogical] = trialSort(obj,directions,NaN,NaN,NaN);
                counts = zeros(1,sum(condLogical));
            end
        end
        
        function counts = dirConditionedCounts(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'win',[176,225])
            addParameter(Parser,'dirs',NaN)
            
            parse(Parser,obj,varargin{:})
            
            obj = Parser.Results.obj;
            win = Parser.Results.win;
            dirs = Parser.Results.dirs;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end       
            
            % Condition on sequences
            counts = nan([1000,length(dirs),length(obj.unitIndex)]);
            maxtrials = 0;
            for di = 1:length(dirs)
                countsTemp = permute(conditionalCounts(obj,win,...
                    dirs(di)),[3,2,1]);
                if size(countsTemp,1) > 1
                    counts(1:size(countsTemp,1),di,:) = countsTemp;
                    maxtrials = max([maxtrials size(countsTemp,1)]);
                else
                    counts(1:length(countsTemp),di,:) = countsTemp;
                    maxtrials = max([maxtrials length(countsTemp)]);
                end
            end
            counts = counts(1:maxtrials,:,:,:);
        end
        
        function sacPrefRaster(obj,dirs,units)
            figure('Name','Direction responses','Position',[559 92 1699 1205]);
            si = 0;
            r = 1;
            xlocs = -r:r;
            ylocs = fliplr(-r:r);
            for xloci = 1:3
                for yloci = 1:3
                    si = si + 1;
                    h(si) = subplot(3,3,si);
                    set(h(si),'visible','off')
                    xy(si,:) = [ylocs(xloci),xlocs(yloci)];
                end
            end
            for di = 1:length(dirs)
                [condInds,condLogical] = trialSort(obj,dirs(di),NaN,NaN,NaN);
                xloc = r*round(cosd(dirs(di)));
                yloc = r*round(sind(dirs(di)));
                si = find(xy(:,1) == yloc & xy(:,2) == xloc);
                set(h(si),'visible','on')
                subplot(h(si))
                t_offsets = [obj.eye(condLogical).targetAcquired];
                rasterPlot(obj,condLogical,units,t_offsets)
                axis tight
                
                ax(di,:) = axis;
                trialN(di) = sum(condLogical);
                count(di) = 0;
                for triali = 1:length(condInds)
                    spikeTimes = obj.spikeTimes{condInds(triali)}{1} - t_offsets(triali);
                    spikeTimes = spikeTimes(ismember(obj.spikeTimes{condInds(triali)}{2},units));
                    count(di) = count(di) + sum(spikeTimes > -300 & spikeTimes < 100);
                end
            end
            lineProps.color = [1 0 0];
            for di = 1:length(dirs)
                xloc = r*round(cosd(dirs(di)));
                yloc = r*round(sind(dirs(di)));
                si = find(xy(:,1) == yloc & xy(:,2) == xloc);
                subplot(h(si))
                axis([min(ax(:,1)) max(ax(:,2)) min(0) max(trialN)])
                plotVertical(0,'lineProperties',lineProps);
                xlabel('Time target onset (ms)')
                ylabel('Trial #')
                mymakeaxis(gca,'xytitle',num2str(dirs(di)))
            end
            si = find(xy(:,1) == 0 & xy(:,2) == 0);
            subplot(h(si))
            polar((pi/180)*([dirs dirs(1)]),([count count(1)])./([trialN trialN(1)]))
        end
        
    end
end