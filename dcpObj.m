%% dcpObj
%
%   Defines properties and methods of object used for analysis of
%   DynamicCoherencePhysiology data.
%
%%

classdef dcpObj
    % DynamicCoherencePhysiology offline analysis class
    properties
        sname;
        datapath;
        spikesExtracted = false;
        trials = 0;
        calib
        unitIndex
        
        dirPref;
        
        speedPref;
        
        initiateCoh;
        
        washout;
        
        dynamicCoh;
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
        
        %% dirPrefTrials
        function obj = dirPrefTrials(obj,trialNs)
        % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:4),obj.sname)
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
                            obj.dirPref.trailtype(ind,1) = ...
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
                            
                            % Add spike times
                            obj.dirPref.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);                            
                            
                        end
                        
                    end
                end
        end
        
        
        %% speedPrefTrials
        function obj = speedPrefTrials(obj,trialNs)
        % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:4),obj.sname)
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
                            obj.speedPref.trialNumbers(ind,1) = ti;
                            obj.speedPref.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.speedPref.trailtype(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.speedPref.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'s\d{3}');
                            obj.speedPref.speeds(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'x\d{3}');
                            obj.speedPref.locations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            [startIndex,endIndex] = regexp(trialname,'y\d{3}');
                            obj.speedPref.locations(ind,2) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            % Add eye information
                            obj.speedPref.eye(:,ind).hpos = (file.data(1,:) - ...
                                mean(file.data(1,obj.calib.t)))*obj.calib.posGain;
                            obj.speedPref.eye(:,ind).vpos = (file.data(2,:) - ...
                                mean(file.data(2,obj.calib.t)))*obj.calib.posGain;
                            
                            obj.speedPref.eye(:,ind).hvel = (file.data(3,:) - ...
                                mean(file.data(3,obj.calib.t)))*obj.calib.speedGain;
                            obj.speedPref.eye(:,ind).vvel = (file.data(4,:) - ...
                                mean(file.data(4,obj.calib.t)))*obj.calib.speedGain;
                            
                            % Add spike times
                            obj.speedPref.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);                            
                            
                        end
                        
                    end
                end
        end
        
        
        %% initiateCohTrials
        function obj = initiateCohTrials(obj,trialNs)
        % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:4),obj.sname)
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
                            % Update trial
                            obj.initiateCoh.trialNumbers(ind,1) = ti;
                            obj.initiateCoh.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.initiateCoh.trailtype(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'h\d{3}');
                            obj.initiateCoh.coh(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.initiateCoh.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'s\d{3}');
                            obj.initiateCoh.speeds(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'x\d{3}');
                            obj.initiateCoh.locations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            [startIndex,endIndex] = regexp(trialname,'y\d{3}');
                            obj.initiateCoh.locations(ind,2) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            % Add eye information
                            obj.initiateCoh.eye(:,ind).hpos = (file.data(1,:) - ...
                                mean(file.data(1,obj.calib.t)))*obj.calib.posGain;
                            obj.initiateCoh.eye(:,ind).vpos = (file.data(2,:) - ...
                                mean(file.data(2,obj.calib.t)))*obj.calib.posGain;
                            
                            obj.initiateCoh.eye(:,ind).hvel = (file.data(3,:) - ...
                                mean(file.data(3,obj.calib.t)))*obj.calib.speedGain;
                            obj.initiateCoh.eye(:,ind).vvel = (file.data(4,:) - ...
                                mean(file.data(4,obj.calib.t)))*obj.calib.speedGain;
                            
                            % Add spike times
                            obj.initiateCoh.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);                            
                            
                        end
                        
                    end
                end
        end
        
        %% dynamicCohTrials
        function obj = dynamicCohTrials(obj,trialNs)
        % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:4),obj.sname)
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
                            obj.dynamicCoh.trialNumbers(ind,1) = ti;
                            obj.dynamicCoh.trialDataFiles{ind} = files(ti+fileInx-1).name;
                            
                            % Parse trial info
                            [startIndex,endIndex] = regexp(trialname,'t\d{3}');
                            obj.dynamicCoh.trailtype(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'q\d{3}');
                            obj.dynamicCoh.sequences(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'p\d{3}');
                            obj.dynamicCoh.perturbations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'d\d{3}');
                            obj.dynamicCoh.directions(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'s\d{3}');
                            obj.dynamicCoh.speeds(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex));
                            
                            [startIndex,endIndex] = regexp(trialname,'x\d{3}');
                            obj.dynamicCoh.locations(ind,1) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            [startIndex,endIndex] = regexp(trialname,'y\d{3}');
                            obj.dynamicCoh.locations(ind,2) = ...
                                str2double(trialname(startIndex+1:endIndex))-100;
                            
                            % Add eye information
                            obj.dynamicCoh.eye(:,ind).hpos = (file.data(1,:) - ...
                                mean(file.data(1,obj.calib.t)))*obj.calib.posGain;
                            obj.dynamicCoh.eye(:,ind).vpos = (file.data(2,:) - ...
                                mean(file.data(2,obj.calib.t)))*obj.calib.posGain;
                            
                            obj.dynamicCoh.eye(:,ind).hvel = (file.data(3,:) - ...
                                mean(file.data(3,obj.calib.t)))*obj.calib.speedGain;
                            obj.dynamicCoh.eye(:,ind).vvel = (file.data(4,:) - ...
                                mean(file.data(4,obj.calib.t)))*obj.calib.speedGain;
                            
                            % Add spike times
                            obj.dynamicCoh.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);                            
                            
                        end
                        
                    end
                end
        end
        
        %% washoutTrials
        function obj = washoutTrials(obj,trialNs)
        % Add direction preference trials to data object
            files = dir(obj.datapath);
            
            % Determine the index of the first data file
            for fileInx = 1:length(files)
                if length(files(fileInx).name) >= length(obj.sname) && ...
                        strcmp(files(fileInx).name(1:4),obj.sname)
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
                            obj.washout.trailtype(ind,1) = ...
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
                            
                            % Add spike times
                            obj.washout.spikeTimes{ind} = ...
                                file.sortedSpikes(obj.unitIndex);                            
                            
                        end
                        
                    end
                end
        end
        
        %% Analysis methods
        function [condInds, condLogical] = dirPrefSort(obj,directions,speeds,locations)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if isnan(directions)
                dMask = true(size(obj.dirPref.directions));
            else
                dMask = ismember(obj.dirPref.directions,directions);
            end
            
            if isnan(speeds)
                sMask = true(size(obj.dirPref.speeds));
            else
                sMask = ismember(obj.dirPref.speeds,speeds);
            end
            
            if any(isnan(locations))
                lMask = true(size(obj.dirPref.locations));
            else
                for li = 1:size(locations,2)
                    lMask(:,li) = ismember(...
                        obj.dirPref.locations(:,li),locations(:,li));
                end
            end
            
            condLogical = prod([dMask,sMask,lMask]);
            condInds = find(condLogical);
        end
        
        
        function [condInds, condLogical] = speedPrefSort(obj,directions,speeds,locations)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if isnan(directions)
                dMask = true(size(obj.speedPref.directions));
            else
                dMask = ismember(obj.speedPref.directions,directions);
            end
            
            if isnan(speeds)
                sMask = true(size(obj.speedPref.speeds));
            else
                sMask = ismember(obj.speedPref.speeds,speeds);
            end
            
            if isnan(cohs)
                sMask = true(size(obj.speedPref.speeds));
            else
                sMask = ismember(obj.speedPref.speeds,speeds);
            end
            
            if any(isnan(locations))
                lMask = true(size(obj.speedPref.locations));
            else
                for li = 1:size(locations,2)
                    lMask(:,li) = ismember(...
                        obj.speedPref.locations(:,li),locations(:,li));
                end
            end
            
            condLogical = prod([dMask,sMask,lMask]);
            condInds = find(condLogical);
        end
        
        
        function [condInds, condLogical] = initiateCohSort(obj,directions,speeds,locations,cohs)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if isnan(directions)
                dMask = true(size(obj.initiateCoh.directions));
            else
                dMask = ismember(obj.initiateCoh.directions,directions);
            end
            
            if isnan(speeds)
                sMask = true(size(obj.initiateCoh.speeds));
            else
                sMask = ismember(obj.initiateCoh.speeds,speeds);
            end
            
            if isnan(cohs)
                cMask = true(size(obj.initiateCoh.coh));
            else
                cMask = ismember(obj.initiateCoh.coh,cohs);
            end
            
            if any(isnan(locations))
                lMask = true(size(obj.initiateCoh.locations));
            else
                lMask = true(size(obj.initiateCoh.locations));
                for li = 1:size(locations,2)
                    lMask(:,li) = ismember(...
                        obj.initiateCoh.locations(:,li),locations(:,li));
                end
            end
            
            condLogical = prod([dMask,sMask,cMask,lMask]);
            condInds = find(condLogical);
        end
        
        function [condInds, condLogical] = dynamicCohSort(obj,directions,speeds,locations,seqs,perts)
        % Find indices of all trials with direction  in directions, speed
        % in speeds, location in locations.
            if isnan(directions)
                dMask = true(size(obj.dynamicCoh.directions));
            else
                dMask = ismember(obj.dynamicCoh.directions,directions);
            end
            
            if isnan(speeds)
                sMask = true(size(obj.dynamicCoh.speeds));
            else
                sMask = ismember(obj.dynamicCoh.speeds,speeds);
            end
            
            if isnan(seqs)
                cMask = true(size(obj.dynamicCoh.coh));
            else
                cMask = ismember(obj.dynamicCoh.coh,cohs);
            end
            
            if isnan(perts)
                pMask = true(size(obj.dynamicCoh.coh));
            else
                pMask = ismember(obj.dynamicCoh.coh,cohs);
            end
            
            if any(isnan(locations))
                lMask = true(size(obj.dynamicCoh.locations));
            else
                lMask = true(size(obj.dynamicCoh.locations));
                for li = 1:size(locations,2)
                    lMask(:,li) = ismember(...
                        obj.dynamicCoh.locations(:,li),locations(:,li));
                end
            end
            
            condLogical = prod([dMask,sMask,cMask,pMask,lMask]);
            condInds = find(condLogical);
        end
                
            
        
        %% Plotting methods
        function h = rasterPlot(obj,trials,units)
        % Raster plot
            h = figure;
            lineProps.Color = 'k';
            lineProps.LineStyle = '-';
            for triali = trials
                plotVertical(obj.dirPref.spikeTimes{triali}{units},...
                    'MinMax',[triali,triali+1],'lineProperties',lineProps);
                hold on
            end
        end
        
    end
    
end