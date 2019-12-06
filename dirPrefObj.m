%% dirPrefObj
%
%   Defines properties and methods of direction preference object used for 
%   analysis of DynamicCoherencePhysiology data.
%
%   Defined as a subclass of the dcpObj
%
%%

classdef dirPrefObj < dcpObj
    properties
        directionTuning;
    end
    
    methods
        %% Core methods
        function obj = dirPrefObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
            obj.directionTuning.tuningFun = @vonMises;
            obj.directionTuning.parameters = nan(1,3);
        end      
        
        
        %% dirPrefTrials
        function obj = dirPrefTrials(obj,trialNs)
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
                        if strcmp(trialname(1:5),'dPref')
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
        
        function dirPrefRaster(obj,dirs,units)
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
                rasterPlot(obj,condLogical,units)
                axis tight
                
                ax(di,:) = axis;
                trialN(di) = sum(condLogical);
                count(di) = 0;
                for triali = 1:length(condInds)
                    spikeTimes = obj.spikeTimes{condInds(triali)}{1};
                    spikeTimes = spikeTimes(ismember(obj.spikeTimes{condInds(triali)}{2},units));
                    count(di) = count(di) + numel(spikeTimes);
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
                xlabel('Time since motion onset (ms)')
                ylabel('Trial #')
                mymakeaxis(gca,'xytitle',num2str(dirs(di)))
            end
            si = find(xy(:,1) == 0 & xy(:,2) == 0);
            subplot(h(si))
            polar((pi/180)*([dirs dirs(1)]),([count count(1)])./([trialN trialN(1)]))
        end
    end
end