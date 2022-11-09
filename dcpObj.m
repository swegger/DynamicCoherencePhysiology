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
        chansIndex
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
    end
     
    properties (SetAccess = private)
        probe;
        unitTypes;
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
        
        function obj = extractKKData(obj,kkname,kkdir,plxname,plxdir,maestrodir,addSpike,startChannel,acceptMU,checkReceipt)
        % Add spiking data to Maestro data files
            if ~exist('addSpike','var')
                addSpike = false;
            end
            if ~exist('startChannel','var')
                startChannel = 'EVT01';
            end
            if ~exist('acceptMU','var')
                acceptMU = false;
            end
            if ~exist('checkReceipt')
                checckReceipt = false;
            end
            
            if checkReceipt
                receiptFile = dir([maestrodir '/spikesReceipt']);
                if isempty(receiptFile)
                    receiptExists = false;
                else
                    receiptExists = true;
                end
            end
            
            if obj.spikesExtracted
                disp('Spikes already extracted and inserted into Maestro files.')
            elseif checkReceipt && receiptExists
                disp('Spike receipt exists, spikes already extracted and inserted into Maestro files.')
                obj.spikesExtracted = true;
                
                tsvFiles = dir([kkdir '/*.tsv']);
                if length(tsvFiles) > 1
                    tsv = tdfread([kkdir '/cluster_info_' kkname(end-5) '.tsv']);
                else
                    tsv = tdfread([kkdir '/cluster_info.tsv']);
                end
                kkIndex = [];
                for uniti = 1:size(tsv.group,1)
                    unitType = [tsv.group(uniti,:) '     '];
                    if acceptMU
                        if strcmp(unitType(1:4),'good') | strcmp(unitType(1:3),'mua')
                            kkIndex = [kkIndex; tsv.channel(uniti) tsv.id(uniti)];
                        end
                    else
                        if strcmp(unitType(1:4),'good')
                            kkIndex = [kkIndex; tsv.channel(uniti) tsv.id(uniti)];
                        end
                    end
                end
                obj.klustaID = kkIndex;
%                 obj = unitsIndex(obj);
%                 obj = setUnitType(obj,kkdir,kkname,acceptMU);
            else
                if strcmp(plxname(end-2:end),'pl2')
                    [~, ~, unitsData, ~, kkIndex] = ...
                        extractKKSpikes(kkname,kkdir,plxname,plxdir,maestrodir,'addSpike',addSpike,...
                        'startChannel',startChannel,'acceptMU',acceptMU);
                    obj.klustaID = kkIndex;
%                     obj = unitsIndex(obj);                    
%                     obj = setUnitType(obj,kkdir,kkname,acceptMU);
                    
                    % Write receipt
                    if ~isempty(kkIndex)
                        Ntrials = length(unitsData);
                        Nspikes = numel(vertcat(unitsData{:,1}));
                        Nclusters = length(kkIndex);
                        Nchannels = numel(unique(vertcat(unitsData{:,3})));
                        ftemp = fopen([maestrodir '/spikesReceipt'],'w+');
                        fprintf(ftemp,'%6s %12s %18s %24s\n','Ntrials','Nspikes','Nclusters','Nchannels');
                        fprintf(ftemp,'%6d %12d %18d %24d\n',Ntrials,Nspikes,Nclusters,Nchannels);
                        fclose(ftemp);
                    else
                        Ntrials = length(unitsData);
                        ftemp = fopen([maestrodir '/spikesReceipt'],'w+');
                        fprintf(ftemp,'%6s %12s %18s %24s\n','Ntrials','Nspikes','Nclusters','Nchannels');
                        fprintf(ftemp,'%6d %12d %18d %24d\n',Ntrials,0,0,0);
                        fclose(ftemp);
                    end
                else
                    error(['Extraction format ' plxname(end-2:end) ' not recognized!'])
                end
                obj.spikesExtracted = true;
            end
        end
        
        function obj = setUnitType(obj,kkdir,kkname)
            if obj.spikesExtracted
                
                if isempty(obj.unitIndex)
                    obj = unitsIndex(obj);
                end
                
                tsvFiles = dir([kkdir '/*.tsv']);
                if length(tsvFiles) > 1
                    tsv = tdfread([kkdir '/cluster_info_' kkname(end-5) '.tsv']);
                else
                    tsv = tdfread([kkdir '/cluster_info.tsv']);
                end
                for uniti = 1:length(obj.unitIndex)
                    unitj = find(obj.unitIndex(uniti) == tsv.id);
                    obj.unitTypes{uniti} = strtrim(tsv.group(unitj,:));
                end
            else
                obj.unitTypes = {''};
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
                chansList = [];
                for ti = fileInx:length(files)
                    file = readcxdata([obj.datapath '/' files(ti).name]);
                    if iscell(file.sortedSpikes)
                        indsList = [indsList file.sortedSpikes{2}];
                        chansList = [chansList file.sortedSpikes{3}];
%                         indsList = [indsList find(~cellfun(@isempty,file.sortedSpikes))];
                    end
                end
                [obj.unitIndex, tempIndx] = unique(indsList);
                obj.chansIndex = chansList(tempIndx);
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
            addParameter(Parser,'t_offsets',NaN)
            
            parse(Parser,obj,width,varargin{:})
            
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            units = Parser.Results.units;
            t = Parser.Results.t;
            trialN = Parser.Results.trialN;     
            t_offsets = Parser.Results.t_offsets;
            
            if any(isnan(trialN))
                trialN = 1:numel(obj.spikeTimes);
            end
            if any(isnan(units))
                units = obj.unitIndex;
%                 if ~isempty(obj.klustaID)
%                     units = obj.klustaID;
%                 else
%                     units = [];
%                     for triali = trialN
%                         units = [units obj.spikeTimes{triali}{2}];
%                     end
%                     units = unique(units);
% %                     units = 1:numel(obj.unitIndex);
%                 end
            end
            t = t(:); 
            
            if isnan(t_offsets)
                t_offsets = zeros(length(trialN),1);
            end
            
            % Find smoothed rates
            f = @(diffs)obj.myBoxCar(diffs,width);
            r = nan(length(t),numel(trialN),numel(units));
            for triali = trialN
                spikeTimes = obj.spikeTimes{triali}(1:2);
                spikeTimes{1} = spikeTimes{1}(ismember(spikeTimes{2},units))-t_offsets(triali);
                spikeTimes{2} = spikeTimes{2}(ismember(spikeTimes{2},units));
                [~, ~, ~, rTemp] = spikeTimes2Rate(spikeTimes,...
                    'time',t,'resolution',1,'Filter',f,...
                    'ComputeVariance','Yes','mixedUnits',true,'units',units);
                r(:,triali,:) = rTemp/(2*width);
%                 [~, ~, ~, rTemp] = spikeTimes2Rate(obj.spikeTimes{triali}(units),...
%                     'time',t,'resolution',1,'Filter',f,...
%                     'ComputeVariance','Yes');
%                 r(:,triali,:) = rTemp/width;
            end
        end
        
        function counts = countSpikes(obj,win,varargin)
        % Calculate trial by trial counts in time window win
            
            % Parse inputs
            Parser = inputParser;
            
            addRequired(Parser,'obj')
            addRequired(Parser,'win')
            addParameter(Parser,'units',NaN)
            addParameter(Parser,'t',-100:1600)
            addParameter(Parser,'trialN',NaN)
            
            parse(Parser,obj,win,varargin{:})
            
            obj = Parser.Results.obj;
            win = Parser.Results.win;
            units = Parser.Results.units;
            t = Parser.Results.t;
            trialN = Parser.Results.trialN;            
            
            if any(isnan(trialN))
                trialN = 1:numel(obj.spikeTimes);
            end
            if any(isnan(units))
                units = obj.unitIndex;
%                 if ~isempty(obj.klustaID)
%                     units = obj.klustaID;
%                 else
%                     units = [];
%                     for triali = trialN
%                         units = [units obj.spikeTimes{triali}{2}];
%                     end
%                     units = unique(units);
% %                     units = 1:numel(obj.unitIndex);
%                 end
            end 
            
            % Find spikeCounts
            for triali = trialN
                spikeTimes = obj.spikeTimes{triali}(1:2);
                spikeTimes{1} = spikeTimes{1}(ismember(spikeTimes{2},units));
                spikeTimes{2} = spikeTimes{2}(ismember(spikeTimes{2},units));
                for uniti = 1:length(units)
                    timeTemp = spikeTimes{1}(spikeTimes{2} == units(uniti));
                    counts(uniti,triali) = sum(timeTemp >= win(1) & timeTemp <= win(2));
                end
            end
        end
        
        function counts = conditionalCounts(obj,win,directions,speeds)
            countsAll = obj.countSpikes(win);
            if exist('countsAll','var')
                [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                    NaN,NaN);
                counts = countsAll(:,condLogical);
            else
                [~,condLogical] = trialSort(obj,directions,speeds,NaN,...
                    NaN,NaN);
                counts = zeros(1,sum(condLogical));
            end
        end
        
        function [r,rste] = conditionalRates(obj,width,directions,speeds,t)
            if ~exist('t','var')
                t = -100:1600;
            end
            rAll = obj.calcRates(width,'t',t);
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                NaN,NaN);
            r = mean(rAll(:,condLogical,:),2);
            rste = sqrt(var(rAll(:,condLogical,:),[],2)/sum(condLogical));
        end
        
        function obj = tableImport(obj)
            fileName = [obj.sname obj.datapath(end-6:end)];
            tableFile = [obj.datapath(1:end-9) ...
                'tables/' obj.sname obj.datapath(end-6:end-1) '.xlsx'];
            T = readtable(tableFile,'Format','auto');
            if str2num(fileName(3:end-1)) > 200101 && str2num(fileName(3:end-1)) < 210000
                % Determine number of electrodes/channels
                if iscell(T.x1_1)
                    if isempty(T.x1_1{1})
                        liveContact(1,1) = false;
                    else
                        liveContact(1,1) = true;
                    end
                else
                    liveContact(1,1) = false;
                end
                if iscell(T.x1_2)
                    if isempty(T.x1_2{1})
                        liveContact(1,2) = false;
                    else
                        liveContact(1,2) = true;
                    end
                else
                    liveContact(1,2) = false;
                end
                if iscell(T.x1_3)
                    if isempty(T.x1_3{1})
                        liveContact(1,3) = false;
                    else
                        liveContact(1,3) = true;
                    end
                else
                    liveContact(1,3) = false;
                end
                if iscell(T.x1_4)
                    if isempty(T.x1_4{1})
                        liveContact(1,4) = false;
                    else
                        liveContact(1,4) = true;
                    end
                else
                    liveContact(1,4) = false;
                end
                if iscell(T.x2_1)
                    if isempty(T.x2_1{1})
                        liveContact(2,1) = false;
                    else
                        liveContact(2,1) = true;
                    end
                else
                    liveContact(2,1) = false;
                end
                if iscell(T.x2_2)
                    if isempty(T.x2_2{1})
                        liveContact(2,2) = false;
                    else
                        liveContact(2,2) = true;
                    end
                else
                    liveContact(2,2) = false;
                end
                if iscell(T.x2_3)
                    if isempty(T.x2_3{1})
                        liveContact(2,3) = false;
                    else
                        liveContact(2,3) = true;
                    end
                else
                    liveContact(2,3) = false;
                end
                if iscell(T.x2_4)
                    if isempty(T.x2_4{1})
                        liveContact(2,4) = false;
                    else
                        liveContact(2,4) = true;
                    end
                else
                    liveContact(2,4) = false;
                end
                if iscell(T.x3_1)
                    if isempty(T.x3_1{1})
                        liveContact(3,1) = false;
                    else
                        liveContact(3,1) = true;
                    end
                else
                    liveContact(3,1) = false;
                end
                if iscell(T.x3_2)
                    if isempty(T.x3_2{1})
                        liveContact(3,2) = false;
                    else
                        liveContact(3,2) = true;
                    end
                else
                    liveContact(3,2) = false;
                end
                if iscell(T.x3_3)
                    if isempty(T.x3_3{1})
                        liveContact(3,3) = false;
                    else
                        liveContact(3,3) = true;
                    end
                else
                    liveContact(3,3) = false;
                end
                if iscell(T.x3_4)
                    if isempty(T.x3_4{1})
                        liveContact(3,4) = false;
                    else
                        liveContact(3,4) = true;
                    end
                else
                    liveContact(3,4) = false;
                end
                liveTrode = sum(liveContact,2)>1;
                
                delims = regexpi(T.Location_x_y_z_{1},',');
                if liveTrode(1)
                    indx = find(strcmp(fileName,T.Date));
                    obj.location.x(1) = str2num(T.Location_x_y_z_{1}(1:delims(1)-1));
                    obj.location.y(1) = str2num(T.Location_x_y_z_{1}(delims(1)+1:delims(2)-1))-0.305;
                    obj.location.z(1) = str2num(T.Location_x_y_z_{1}(delims(2)+1:end));
                    obj.location.depth(1) = str2num(T.Location_x_y_z_{indx-1})-...
                        str2num(T.Location_x_y_z_{find(strcmp(T.Location_x_y_z_,'Depth 1'))+1});
                else
                    obj.location.x(1) = NaN;
                    obj.location.y(1) = NaN;
                    obj.location.z(1) = NaN;
                    obj.location.depth(1) = NaN;
                end
                
                if liveTrode(2)
                    indx = find(strcmp(fileName,T.x1_2));
                    obj.location.x(2) = str2num(T.Location_x_y_z_{1}(1:delims(1)-1));
                    obj.location.y(2) = str2num(T.Location_x_y_z_{1}(delims(1)+1:delims(2)-1));
                    obj.location.z(2) = str2num(T.Location_x_y_z_{1}(delims(2)+1:end));
                    obj.location.depth(2) = str2num(T.x1_3{indx-1})-...
                        str2num(T.x1_3{find(strcmp(T.x1_3,'Depth 2'))+1});
                else
                    obj.location.x(2) = NaN;
                    obj.location.y(2) = NaN;
                    obj.location.z(2) = NaN;
                    obj.location.depth(2) = NaN;
                end
                
                if liveTrode(3)
                    indx = find(strcmp(fileName,T.x2_3));
                    obj.location.x(3) = str2num(T.Location_x_y_z_{1}(1:delims(1)-1));
                    obj.location.y(3) = str2num(T.Location_x_y_z_{1}(delims(1)+1:delims(2)-1))+0.305;
                    obj.location.z(3) = str2num(T.Location_x_y_z_{1}(delims(2)+1:end));
                    obj.location.depth(3) = str2num(T.x2_4{indx-1})-...
                        str2num(T.x2_4{find(strcmp(T.x2_4,'Depth 3'))+1});
                else
                    obj.location.x(3) = NaN;
                    obj.location.y(3) = NaN;
                    obj.location.z(3) = NaN;
                    obj.location.depth(3) = NaN;
                end
                
            elseif str2num(fileName(3:end-1)) >= 210000
                if strcmp(T.Probe{1},'24V')
                    chanMap = [7 0; 6 1; 5 2; 4 3; 3 4; 2 5; 1 6; 0 7; 23 8; 22 9; 21 10; 20 11; 19 12; 18 13; 17 14; 16 15; 15 16; 14 17; 13 18; 12 19; 11 20; 10 21; 9 22; 8 23];
                    indx = find(strcmp(fileName,T.Date));
                    tempLoc = regexp(T.Location_x_y_z_{1},',','split');
                    obj.location.x = repmat(str2num(tempLoc{1}),[1,24]);
                    obj.location.y = repmat(str2num(tempLoc{2}),[1,24]);
                    obj.location.z = repmat(str2num(tempLoc{3}),[1,24]);
%                     obj.location.x = repmat(str2num(T.Location_x_y_z_{1}(1:4)),[1,24]);
%                     obj.location.y = repmat(str2num(T.Location_x_y_z_{1}(6:9)),[1,24]);
%                     obj.location.z = repmat(str2num(T.Location_x_y_z_{1}(11:14)),[1,24]);
                    obj.location.depth = nan(1,24);
                    refDepth = str2num(T.Location_x_y_z_{indx-1})-...
                        str2num(T.Location_x_y_z_{find(strcmp(T.Location_x_y_z_,'Depth'))+1});
                    chans = unique(obj.chansIndex);
                    for chani = 1:length(chans)
                        chanLoc = chanMap(chanMap(:,1) == chans(chani),2);
                        obj.location.depth(chanLoc+1) = refDepth-150*chanLoc;
                    end
                elseif strcmp(T.Probe{1},'24V.2')
                    chanMap = [23 0; 22 1; 21 2; 20 3; 19 4; 18 5; 17 6; 16 7; 15 8; 14 9; 13 10; 12 11; 11 12; 10 13; 9 14; 8 15; 7 16; 6 17; 5 18; 4 19; 3 20; 2 21; 1 22; 0 23];
                    indx = find(strcmp(fileName,T.Date));
                    tempLoc = regexp(T.Location_x_y_z_{1},',','split');
                    obj.location.x = repmat(str2num(tempLoc{1}),[1,24]);
                    obj.location.y = repmat(str2num(tempLoc{2}),[1,24]);
                    obj.location.z = repmat(str2num(tempLoc{3}),[1,24]);
%                     obj.location.x = repmat(str2num(T.Location_x_y_z_{1}(1:4)),[1,24]);
%                     obj.location.y = repmat(str2num(T.Location_x_y_z_{1}(6:9)),[1,24]);
%                     obj.location.z = repmat(str2num(T.Location_x_y_z_{1}(11:14)),[1,24]);
                    obj.location.depth = nan(1,24);
                    refDepth = str2num(T.Location_x_y_z_{indx-1})-...
                        str2num(T.Location_x_y_z_{find(strcmp(T.Location_x_y_z_,'Depth'))+1});
                    chans = unique(obj.chansIndex);
                    for chani = 1:length(chans)
                        chanLoc = chanMap(chanMap(:,1) == chans(chani),2);
                        obj.location.depth(chanLoc+1) = refDepth-150*chanLoc;
                    end
                elseif strcmp(T.Probe{1},'single')
                    chanMap = [0,0];
                    indx = find(strcmp(fileName,T.Date));
                    tempLoc = regexp(T.Location_x_y_z_{1},',','split');
                    obj.location.x = str2num(tempLoc{1});
                    obj.location.y = str2num(tempLoc{2});
                    obj.location.z = str2num(tempLoc{3});
                    refDepth = str2num(T.Location_x_y_z_{indx-1})-...
                        str2num(T.Location_x_y_z_{find(strcmp(T.Location_x_y_z_,'Depth'))+1});
                   obj.location.depth = refDepth;
                end
            else
                indx = find(strcmp(fileName,T.Date));
                obj.location.x = str2num(T.Location_x_y_z_{1}(1:4));
                obj.location.y = str2num(T.Location_x_y_z_{1}(6:9));
                obj.location.z = str2num(T.Location_x_y_z_{1}(11:14));
                obj.location.depth = str2num(T.Location_x_y_z_{indx-1})-...
                    str2num(T.Location_x_y_z_{find(strcmp(T.Location_x_y_z_,'Depth'))+1});
            end
            
        end
        
        function obj = addProbeInfo(obj)
            fileName = [obj.sname obj.datapath(end-6:end)];
            tableFile = [obj.datapath(1:end-9) ...
                'tables/' obj.sname obj.datapath(end-6:end-1) '.xlsx'];
            T = readtable(tableFile,'Format','auto');
            if str2num(fileName(3:end-1)) > 200101 && str2num(fileName(3:end-1)) < 210000
                % Determine number of electrodes/channels
                if iscell(T.x1_1)
                    if isempty(T.x1_1{1})
                        liveContact(1,1) = false;
                    else
                        liveContact(1,1) = true;
                    end
                else
                    liveContact(1,1) = false;
                end
                if iscell(T.x1_2)
                    if isempty(T.x1_2{1})
                        liveContact(1,2) = false;
                    else
                        liveContact(1,2) = true;
                    end
                else
                    liveContact(1,2) = false;
                end
                if iscell(T.x1_3)
                    if isempty(T.x1_3{1})
                        liveContact(1,3) = false;
                    else
                        liveContact(1,3) = true;
                    end
                else
                    liveContact(1,3) = false;
                end
                if iscell(T.x1_4)
                    if isempty(T.x1_4{1})
                        liveContact(1,4) = false;
                    else
                        liveContact(1,4) = true;
                    end
                else
                    liveContact(1,4) = false;
                end
                if iscell(T.x2_1)
                    if isempty(T.x2_1{1})
                        liveContact(2,1) = false;
                    else
                        liveContact(2,1) = true;
                    end
                else
                    liveContact(2,1) = false;
                end
                if iscell(T.x2_2)
                    if isempty(T.x2_2{1})
                        liveContact(2,2) = false;
                    else
                        liveContact(2,2) = true;
                    end
                else
                    liveContact(2,2) = false;
                end
                if iscell(T.x2_3)
                    if isempty(T.x2_3{1})
                        liveContact(2,3) = false;
                    else
                        liveContact(2,3) = true;
                    end
                else
                    liveContact(2,3) = false;
                end
                if iscell(T.x2_4)
                    if isempty(T.x2_4{1})
                        liveContact(2,4) = false;
                    else
                        liveContact(2,4) = true;
                    end
                else
                    liveContact(2,4) = false;
                end
                if iscell(T.x3_1)
                    if isempty(T.x3_1{1})
                        liveContact(3,1) = false;
                    else
                        liveContact(3,1) = true;
                    end
                else
                    liveContact(3,1) = false;
                end
                if iscell(T.x3_2)
                    if isempty(T.x3_2{1})
                        liveContact(3,2) = false;
                    else
                        liveContact(3,2) = true;
                    end
                else
                    liveContact(3,2) = false;
                end
                if iscell(T.x3_3)
                    if isempty(T.x3_3{1})
                        liveContact(3,3) = false;
                    else
                        liveContact(3,3) = true;
                    end
                else
                    liveContact(3,3) = false;
                end
                if iscell(T.x3_4)
                    if isempty(T.x3_4{1})
                        liveContact(3,4) = false;
                    else
                        liveContact(3,4) = true;
                    end
                else
                    liveContact(3,4) = false;
                end
                obj.probe.type = 'tet';
                obj.probe.liveContacts = liveContact;
                
            elseif str2num(fileName(3:end-1)) >= 210000
                if strcmp(T.Probe{1},'24V')
                    obj.probe.type = '24V';
                    obj.probe.liveContacts = true(1,24);
                end
            else
                obj.probe.type = 'tet';
                obj.probe.liveContacts = [true(1,4); false(3,4)];
            end
        end
        
        function GetSize(this)
            props = properties(this);
            totSize = 0;
            
            for ii=1:length(props)
                currentProperty = getfield(this, char(props(ii)));
                s = whos('currentProperty');
                totSize = totSize + s.bytes;
            end
            
            fprintf(1, '%d bytes\n', totSize);
        end
        
        function mysaveobj(obj,destination)
            s = struct(obj);
            save(destination,'-struct','s')
        end
        
        function [cc, shufflecc] = correlograms(obj,samplerate,unitsIndex,win,shuffleN)
        % Compute differences in spikes times from all recorded spikes
            if ~exist('shuffleN','var')
                shuffleN = NaN;
            end
            dataShort = obj.datapath(end-6:end);
            if strcmp(obj.sname,'ar')
                datapath2 = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' obj.datapath(end-8:end)];
            end
            if str2num(dataShort(1:end-1)) > 191114
                kkfile = [obj.sname dataShort '.kwik'];
                kkpath = [datapath2(1:end-1) 'kk'];
                kkfull = [kkpath '/' kkfile];
                kkflg = exist(kkfull,'file');
                plxflg = false;
            else
                kkflg = false;
                plxfile = [obj.sname dataShort '.plx'];
                plxpath = [datapath2(1:end-1) 'plx'];
                plxfull = [plxpath '/' plxfile];
                plxflg = exist(plxfull,'file');
            end
            
            % Make edges from bin centers specified by win
            edges = [win-(win(2)-win(1))/2 win(end)+(win(2)-win(1))/2];
            
            % Preallocate correlograms
            if isnan(shuffleN)
                shufflecc = nan(length(win),length(unitsIndex),length(unitsIndex),1);
            else
                shufflecc = zeros([length(win),length(unitsIndex),length(unitsIndex),shuffleN]);
            end
            cc = zeros([length(win),length(unitsIndex),length(unitsIndex)]);
            
            if kkflg
                % Get spike times and cell ID
                spktimes = double(hdf5read(kkfull, '/channel_groups/0/spikes/time_samples'))/samplerate;
                clusters = hdf5read(kkfull, '/channel_groups/0/spikes/clusters/main');
                
                spktimes = spktimes(ismember(clusters,unitsIndex));
                clusters = clusters(ismember(clusters,unitsIndex));
            elseif plxflg
                unitNums = size(unitsIndex, 2);
                plx = readPLXFileC(plxfull,'waves');
                for i = 1:unitNums
                    unitsDataTmp{i} = plx.SpikeChannels(obj.chansIndex(obj.unitIndex==unitsIndex(i))).Timestamps(...
                        plx.SpikeChannels(obj.chansIndex(obj.unitIndex==unitsIndex(i))).Units == unitsIndex(i));
                    clustersTmp{i} = unitsIndex(i)*ones(size(unitsDataTmp{i}));
                end
                spktimes = double(vertcat(unitsDataTmp{:}))/samplerate;
                clusters = vertcat(clustersTmp{:});
            end
            
            % Find differences and add to correlogram
            if plxflg || kkflg
                for uniti = 1:length(unitsIndex)
                    for unitj = 1:length(unitsIndex)
                        ts1 = spktimes(clusters == unitsIndex(uniti));
                        ts2 = spktimes(clusters == unitsIndex(unitj));
                        ind = 1;
                        for ti = 1:length(ts1)
                            d = ts1(ti) - ts2;
                            tdiffs = histc(d,edges);
                            cc(:,uniti,unitj) = cc(:,uniti,unitj) + tdiffs(1:end-1);
                            
                            if isnan(shuffleN)
                                shufflediffs(:,uniti,unitj,:) = NaN;
                            else
                                nspikes = sum(d>=edges(1) & d<edges(end));
                                if nspikes == 1
                                    for ri = 1:shuffleN
                                        randtimes = (win(end)-win(1))*rand + win(1);
                                        shufflediffs = histc(randtimes,edges);
                                        shufflecc(:,uniti,unitj,ri) = shufflecc(:,uniti,unitj,ri) + shufflediffs(1:end-1)';
                                    end
                                else
                                    randtimes = (win(end)-win(1))*rand(nspikes,shuffleN) + win(1);
                                    shufflediffs = histc(randtimes,edges); % Preserves rates, eliminates timing
                                    shufflecc(:,uniti,unitj,:) = shufflecc(:,uniti,unitj,:) + permute(shufflediffs(1:end-1,:),[1,4,3,2]);
                                end
                            end
                        end
                    end
                end
            else
                for i = 1:length(unitsIndex)
                    for j = 1:length(unitsIndex)
                        cc = nan([length(win),length(unitsIndex),length(unitsIndex)]);
                        shufflecc = nan([length(win),length(unitsIndex),length(unitsIndex),shuffleN]);
                    end
                end
            end
        
        end
        
        function counts = spikeCount(obj,unitsIndex)
            dataShort = obj.datapath(end-6:end);
            if strcmp(obj.sname,'ar')
                datapath2 = ['/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/' obj.datapath(end-8:end)];
            end
            if str2num(dataShort(1:end-1)) > 191114
                kkfile = [obj.sname dataShort '.kwik'];
                kkpath = [datapath2(1:end-1) 'kk'];
                kkfull = [kkpath '/' kkfile];
                kkflg = exist(kkfull,'file');
                plxflg = false;
            else
                kkflg = false;
                plxfile = [obj.sname dataShort '.plx'];
                plxpath = [datapath2(1:end-1) 'plx'];
                plxfull = [plxpath '/' plxfile];
                plxflg = exist(plxfull,'file');
            end
            
            if kkflg
                % Get spike times and cell ID
                clusters = hdf5read(kkfull, '/channel_groups/0/spikes/clusters/main');
                clusters = clusters(ismember(clusters,unitsIndex));
            elseif plxflg
                unitNums = size(unitsIndex, 2);
                plx = readPLXFileC(plxfull,'waves');
                for i = 1:unitNums
                    unitsDataTmp{i} = plx.SpikeChannels(obj.chansIndex(obj.unitIndex==unitsIndex(i))).Timestamps(...
                        plx.SpikeChannels(obj.chansIndex(obj.unitIndex==unitsIndex(i))).Units == unitsIndex(i));
                    clustersTmp{i} = unitsIndex(i)*ones(size(unitsDataTmp{i}));
                end
                clusters = vertcat(clustersTmp{:});
            end
            
            if plxflg || kkflg
                for uniti = 1:length(unitsIndex)
                    counts(uniti) = sum(clusters == unitsIndex(uniti));
                end
            else
                counts = nan(length(unitsIndex));
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
        function h = rasterPlot(obj,trials,units,t_offsets)
        % Raster plot
%             h = figure;
            if ~exist('t_offsets','var')
                t_offsets = [];
            end
            if islogical(trials)
                trials = find(trials);
            end
            lineProps.Color = 'k';
            lineProps.LineStyle = '-';
            for triali = 1:length(trials)
                spikeTimes = obj.spikeTimes{trials(triali)}{1};
                spikeTimes = spikeTimes(ismember(obj.spikeTimes{trials(triali)}{2},units));
                if ~isempty(spikeTimes) && isempty(t_offsets)
                    plotVertical(spikeTimes,...
                        'MinMax',[triali,triali+1],'lineProperties',lineProps);
                elseif ~isempty(spikeTimes)
                    plotVertical(spikeTimes-t_offsets(triali),...
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
            addParameter(Parser,'plotPatch',true)
            
            parse(Parser,obj,condLogical,varargin{:})
            
            obj = Parser.Results.obj;
            condLogical = Parser.Results.condLogical;
            h = Parser.Results.h;
            sh = Parser.Results.sh;
            t = Parser.Results.t;
            color = Parser.Results.color;
            normalizer = Parser.Results.normalizer;
            plotPatch = Parser.Results.plotPatch;
            
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
            
            if plotPatch
                patchProps.FaceAlpha = 0.3;
                if ~any(isnan(color))
                    patchProps.FaceColor = color;
                else
                    color = [0 0 0];
                    patchProps.FaceColor = color;
                end
                myPatch(t(:),mE(:),steE(:),'patchProperties',patchProps);
            end
            hold on
            plot(t,mE,'Color',color,'LineWidth',2)
        end
        
        
        
        
    end
    
end

