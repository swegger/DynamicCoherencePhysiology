%% speedPrefObj
%
%   Defines properties and methods of speed preference object used for 
%   analysis of DynamicCoherencePhysiology data.
%
%   Defined as a subclass of the dcpObj
%
%%

classdef speedPrefObj < dcpObj
    properties
        speedTuning;
    end
    
    methods
        %% Core methods
        function obj = speedPrefObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
            obj.speedTuning.tuningFun = @vonMises;
            obj.speedTuning.parameters = nan(1,3);
        end 
        
        %% speedPrefTrials
        function obj = speedPrefTrials(obj,trialNs)
        % Add speed preference trials to data object
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
                        if strcmp(trialname(1:5),'sPref')
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
                                    file.sortedSpikes;  
%                                 obj.spikeTimes{ind} = ...
%                                     file.sortedSpikes(obj.unitIndex);                                
                            else
                                obj.spikeTimes{ind}{1} = file.spikes;                    
                            end
                            
                        end
                        
                    end
                end
        end
        
        function speedPrefRaster(obj,directions,speeds,units)
            figure('Name','Speed responses','Position',[559 1011 1699 500]);
            spi = 0;
            for di = 1:length(directions)
                for si = 1:length(speeds)+1
                    spi = spi + 1;
                    h(spi) = subplot(length(directions),length(speeds)+1,spi);
                    spInd(di,si) = spi;
                end
            end
            for di = 1:length(directions)
                for si = 1:length(speeds)
                    [condInds,condLogical] = trialSort(obj,directions(di),speeds(si),NaN,NaN);
                    subplot(h(spInd(di,si)))
                    rasterPlot(obj,condLogical,units)
                    axis tight
                    
                    ax(si,:) = axis;
                    trialN(si) = sum(condLogical);
                    count(di,si) = 0;
                    for triali = 1:length(condInds)
                        spikeTimes = obj.spikeTimes{condInds(triali)}{1};
                        spikeTimes = spikeTimes(ismember(obj.spikeTimes{condInds(triali)}{2},units));
                        count(di,si) = count(si) + numel(spikeTimes);
                    end
                end
            end
            lineProps.color = [1 0 0];
            for di = 1:length(directions)
                for si = 1:length(speeds)
                    subplot(h(spInd(di,si)))
                    axis([min(ax(:,1)) max(ax(:,2)) min(0) max(trialN)])
                    plotVertical(0,'lineProperties',lineProps);
                    xlabel('Time since motion onset (ms)')
                    ylabel('Trial #')
                    mymakeaxis(gca,'xytitle',num2str(speeds(si)))
                end
            end
            for di = 1:length(directions)
                subplot(h(spInd(di,length(speeds)+1)))
                plot(speeds,count(di,:),'-o')
                xlabel('Speed (deg/s)')
                ylabel('Count (spikes)')
                mymakeaxis(gca,'xytitle',['Dir = ' num2str(directions(di))])
            end
        end
        
    end
end