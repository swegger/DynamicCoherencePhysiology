%% dirPref
%
%   Defines properties and methods of direction preference object used for 
%   analysis of DynamicCoherencePhysiology data.
%
%   Defined as a subclass of the dcpObj
%
%%

classdef dirPref < dcpObj
    properties
        trialNumbers;
        trialDataFiles;
        trialtype;
        directions;
        speeds;
        locations;
        eye;
    end
    
    methods
        %% Class constructor
        function obj = dirPrefObj()
            
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
                            obj.dirPref.trialNumbers(ind,1) = ti;
                            obj.dirPref.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.dirPref.trialtpe(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.dirPref.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'s\d{3}');
                            obj.dirPref.speeds(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'x\d{3}');
                            obj.dirPref.locations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            [startIndex,endIndex] = regexp(trialname,'y\d{3}');
                            obj.dirPref.locations(ind,2) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            % Add eye information
                            obj.dirPref.eye(:,ind).hpos = (file.data(1,:) - ...
                                mean(file.data(1,obj.calib.t)))*obj.calib.posGain;
                            obj.dirPref.eye(:,ind).vpos = (file.data(2,:) - ...
                                mean(file.data(2,obj.calib.t)))*obj.calib.posGain;
                            
                            obj.dirPref.eye(:,ind).hvel = (file.data(3,:) - ...
                                mean(file.data(3,obj.calib.t)))*obj.calib.speedGain;
                            obj.dirPref.eye(:,ind).vvel = (file.data(4,:) - ...
                                mean(file.data(4,obj.calib.t)))*obj.calib.speedGain;
                            
                            sacs = saccadeDetect(file.data(3,:)*obj.calib.speedGain,...
                                file.data(4,:)*obj.calib.speedGain,...
                                'accelerationThreshold',obj.calib.accThres,...
                                'windowSize',40);
                            obj.dirPref.eye(:,ind).hvel(sacs) = NaN;
                            obj.dirPref.eye(:,ind).vvel(sacs) = NaN;
                            obj.dirPref.eye(:,ind).saccades = sacs;
                            
                            % Add spike times
                            obj.dirPref.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);                            
                            
                        end
                        
                    end
                end
        end
    end
end