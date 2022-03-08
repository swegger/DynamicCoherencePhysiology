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
        objType = 'dynamicCohObj';
        sequences;
        perturbations;
        r;
        R;
        Rste;
        coh;
        
        preferredDirection;
        preferredDirectionRelative;
        
        Epoch;
        
        saveLocation;
    end
    
    properties (SetAccess = private)
        filterWidth;
        preferredDirectionWin;
        rateCutoff;
        cutWindow;
        passCutoff;
        pertAmp;
        
        eye_t;
        neuron_t;
    end
    
    methods
        %% Core methods
        function obj = dynamicCohObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
        end 
        
        %% dynamicCohTrials
        function obj = dynamicCohTrials(obj,trialNs,forceRead)
            if ~exist('forceRead','var')
                forceRead = false;
            end
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
                            
                            if isempty(file.trialInfo.perts)
                                obj.pertAmp(ind) = 0;
                            else
                                obj.pertAmp(ind) = file.trialInfo.perts.amp;
                            end
                            
                            % Add spike times
                            if obj.spikesExtracted
                                obj.spikeTimes{ind} = ...
                                    file.sortedSpikes;    
%                                 obj.spikeTimes{ind} = ...
%                                     file.sortedSpikes(obj.unitIndex);                             
                            else
                                obj.spikeTimes{ind}{1} = file.spikes;                    
                            end                           
                            
                            onTemp = [file.targets.on{2:end}];
                            t_temp = min(onTemp):(max(onTemp)-1);
                            
                            
                        elseif forceRead
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
                            
                            if isempty(file.trialInfo.perts)
                                obj.pertAmp(ind) = 0;
                            else
                                obj.pertAmp(ind) = file.trialInfo.perts.amp;
                            end
                            
                            % Add spike times
                            if obj.spikesExtracted
                                obj.spikeTimes{ind} = ...
                                    file.sortedSpikes;    
%                                 obj.spikeTimes{ind} = ...
%                                     file.sortedSpikes(obj.unitIndex);                             
                            else
                                obj.spikeTimes{ind}{1} = file.spikes;                    
                            end                           
                            
                            onTemp = [file.targets.on{2:end}];
                            t_temp = min(onTemp):(max(onTemp)-1);
                        end
                        
                    end
                end
                if exist('t_temp','var')
                    obj.eye_t = t_temp;
                else
                    obj.eye_t = [];
                end
        end
        
        %% addCoh
        function obj = addCoh(obj)
            if isempty(obj.eye_t)
                t = -100:1600;
            else
                t = obj.eye_t;
            end
            obj.coh = nan(size(t,2),5);
            
            % Initial
            obj.coh(t <= 450,:) = 60;
            
            % First block
            obj.coh(t > 450 & t <= 750,1) = 100;
            obj.coh(t > 450 & t <= 750,2) = 20;
            obj.coh(t > 450 & t <= 750,3) = 100;
            obj.coh(t > 450 & t <= 750,4) = 20;
            obj.coh(t > 450 & t <= 750,5) = 60;
            
            % Second block            
            obj.coh(t > 750 & t <= 1050,1) = 20;
            obj.coh(t > 750 & t <= 1050,2) = 100;
            obj.coh(t > 750 & t <= 1050,3) = 100;
            obj.coh(t > 750 & t <= 1050,4) = 20;
            obj.coh(t > 750 & t <= 1050,5) = 60;
            
            % Final
            obj.coh(t > 1050,:) = 60;
            
        end
        
        %% Behavior analysis methods
        function meanEyeSpeeds = dynamicCohMeanEyeBoot(obj,dirs,bootN,sampN)
        % Mean eye speed for each sequence and target speed
            speeds = unique(obj.speeds);
            seqs = unique(obj.sequences);
            for si = 1:length(speeds)
                for seqi = 1:length(seqs)
                    condInds = trialSort(obj,dirs,speeds(si),NaN,NaN,seqs(seqi),NaN);
                    for booti = 1:bootN
                        condIndstemp = randsample(condInds,sampN);
                        meanEyeSpeeds(:,seqi,si,booti) = MeanEyeSpeed(obj,condIndstemp);
                    end
                end
            end
        end
        
        function obj = addBehavioralEpochs(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'dirs',[0,180])
            addParameter(Parser,'seqs',1:5)
            addParameter(Parser,'t',0:(1600-1))
            addParameter(Parser,'EpochT',[0 150; 150 450; 450 750; 750 1050; 1050 1350])
            
            parse(Parser,obj,varargin{:})
            
            obj = Parser.Results.obj;
            dirs = Parser.Results.dirs;
            seqs = Parser.Results.seqs;
            EpochT = Parser.Results.EpochT;
            t = Parser.Results.t;
            
            % Find mean eye speeds for each sequence
            for di = 1:length(dirs)
                for seqi = 1:length(seqs)
                    condInds = trialSort(obj,dirs(di),NaN,NaN,NaN,seqs(seqi),NaN);
                    meanEyeSpeeds(:,seqi,di) = MeanEyeSpeed(obj,condInds);
                end
            end
            
            obj.Epoch.eye = nan(max(diff(EpochT,1,2)),length(seqs),4,length(seqs)+1,length(dirs));
            for epochi = 1:size(EpochT,1)
                ttemp = 0:(EpochT(epochi,2)-EpochT(epochi,1));
                for controli = 1:(length(seqs)+1)
                    for seqi = 1:length(seqs)
                        for di = 1:length(dirs)
                            if controli-1 == 0
                                obj.Epoch.eye(1:length(ttemp),seqi,epochi,controli,di) = ...
                                    meanEyeSpeeds(t>=EpochT(epochi,1) & t<= EpochT(epochi,2),seqi,di);
                                
                                obj.Epoch.p(seqi,epochi,controli,di).eye = ...
                                    pieceWiseExpRiseFitter(ttemp(:),...
                                        obj.Epoch.eye(1:length(ttemp),seqi,epochi,controli,di),...
                                        [0,10,200,50],[Inf,Inf,Inf,200],[-Inf,-Inf,0,0]);
                            elseif controli-1 ~= seqi
                                obj.Epoch.eye(1:length(ttemp),seqi,epochi,controli,di) = ...
                                    meanEyeSpeeds(t>=EpochT(epochi,1) & t<= EpochT(epochi,2),seqi,di) -... 
                                    meanEyeSpeeds(t>=EpochT(epochi,1) & t<= EpochT(epochi,2),controli-1,di);
                                
                                obj.Epoch.p(seqi,epochi,controli,di).eye = ...
                                    pieceWiseExpRiseFitter(ttemp(:),...
                                        obj.Epoch.eye(1:length(ttemp),seqi,epochi,controli,di),...
                                        [0,3,200,50],[Inf,Inf,Inf,200],[-Inf,-Inf,0,0]);
                            else
                                obj.Epoch.p(seqi,epochi,controli,di).eye = nan(1,4);
                            end
                        end
                    end
                end
            end
        end   
        
        function [C, res] = findCovarianceBehavior(obj,sequences,perturbations,dirs,win,interpolationMethod)
            if ~exist('interpolationMethod','var')
                interpolationMethod = 'Linear';
            end
            res = nan(1000,length(obj.eye_t(obj.eye_t>=win(1) & obj.eye_t<= win(2))));
            ind = 1;
            for di = 1:length(dirs)
                for seqi = 1:length(sequences)
                    [~,condLogical] = trialSort(obj,dirs(di),NaN,NaN,NaN,...
                        sequences(seqi),perturbations);
                    eh = vertcat(obj.eye(:).hvel);
                    ev = vertcat(obj.eye(:).vvel);
                    eSpeed = sqrt(eh(condLogical,:).^2 + ...
                        ev(condLogical,:).^2);
                    
                    % interpolate during sacceades
                    switch interpolationMethod
                        
                        case {'linear','Linear'}
                            % Linearly interpolate between points
                            for triali = 1:size(eSpeed,1)
                                Vq = interp1(obj.eye_t(~isnan(eSpeed(triali,:))),...
                                    eSpeed(triali,~isnan(eSpeed(triali,:))),...
                                    obj.eye_t(isnan(eSpeed(triali,:))));
                                eSpeed(triali,isnan(eSpeed(triali,:))) = Vq;
                            end
                            
                        case {'GP','GaussianProcess','gaussianProcess','gaussianprocess'}
                            % TO DO
                            
                    end
                    
                    res(ind:ind+sum(condLogical)-1,:) = ...
                        eSpeed(:,obj.eye_t>=win(1) & obj.eye_t<= win(2)) - ...
                        mean(eSpeed(:,obj.eye_t>=win(1) & obj.eye_t<= win(2)),1);
                    ind = ind+sum(condLogical);
                end
            end
            res = res(1:ind-1,:);
            stdRes = std(res(:));
            res = res(~any(abs(res) > 1.5*stdRes,2),:);
            C = cov(res);
            
        end
        
        %% Neural analysis methods
        function [r,rste] = conditionalRates(obj,width,directions,speeds,...
                sequences,perturbations,t)
            if ~exist('t','var')
                t = -100:1600;
            end
            
            % Determine if any spikies were recorded
            units = [];
            for triali = 1:length(obj.spikeTimes)
                if length(obj.spikeTimes{triali}) > 1
                    units = [units obj.spikeTimes{triali}{2}];
                end
            end
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                sequences,perturbations);
            if isempty(units)
                warning(['No spike times recorded for dynamicCohObj associated with ' obj.datapath(end-8:end)])
                r = nan(length(t),1,length(obj.unitIndex));
                rste = nan(length(t),1,length(obj.unitIndex));
            else
                rAll = obj.calcRates(width,'t',t);
                r = mean(rAll(:,condLogical,:),2);
                rste = sqrt(var(rAll(:,condLogical,:),[],2)/sum(condLogical));
            end
        end
        
        function obj = dynamicCohSeqConditionedRates(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'width',50)
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'seqs',NaN)
            addParameter(Parser,'perts',NaN)
            addParameter(Parser,'t',-100:1600)
            addParameter(Parser,'marginalizeDirection',true)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            dirs = Parser.Results.dirs;
            speeds = Parser.Results.speeds;
            locs = Parser.Results.locs;
            seqs = Parser.Results.seqs;
            perts = Parser.Results.perts;
            t = Parser.Results.t;
            marginalizeDirection = Parser.Results.marginalizeDirection;
            
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
                if marginalizeDirection
                    [R(:,seqi,:),Rste(:,seqi,:)] = conditionalRates(obj,width,...
                        dirs,speeds,seqs(seqi),perts,t);
                else
                    for di = 1:length(dirs)
                        [R(:,seqi,:,di),Rste(:,seqi,:,di)] = conditionalRates(obj,width,...
                            dirs(di),speeds,seqs(seqi),perts,t);                    
                    end
                end
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
        
        function counts = conditionalCounts(obj,win,directions,speeds,...
                sequences,perturbations)
            if length(obj.spikeTimes{1}) < 2
                warning(['No spike times recorded for dynamicCohObj associated with ' obj.datapath(end-8:end)])
            else
                countsAll = obj.countSpikes(win);
            end
            if exist('countsAll','var')
                [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                    sequences,perturbations);
                counts = countsAll(:,condLogical);
            else
                [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                    sequences,perturbations);
                counts = zeros(1,sum(condLogical));
            end
        end
        
        function counts = seqConditionedCounts(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'win',[176,225])
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'seqs',NaN)
            addParameter(Parser,'perts',NaN)
            addParameter(Parser,'locs',NaN)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            win = Parser.Results.win;
            dirs = Parser.Results.dirs;
            speeds = Parser.Results.speeds;
            seqs = Parser.Results.seqs;
            perts = Parser.Results.perts;
            locs = Parser.Results.locs;
            
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
            counts = nan(1000,length(seqs),length(obj.unitIndex));
            maxcounts = 0;
            for seqi = 1:length(seqs)
                countsTemp = permute(conditionalCounts(obj,win,...
                    dirs,speeds,seqs(seqi),perts),[2,3,1]);
                if size(countsTemp,1) > 1
                    counts(1:size(countsTemp,1),seqi,:) = countsTemp;
                    maxcounts = max([maxcounts size(countsTemp,1)]);
                else
                    counts(1:length(countsTemp),seqi,:) = countsTemp;
                    maxcounts = max([maxcounts length(countsTemp)]);
                end
            end
            counts = counts(1:maxcounts,:,:);
        end
        
        function obj = addNeuralEpochs(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'dirs',[0,180])
            addParameter(Parser,'seqs',1:5)
            addParameter(Parser,'t',-100:1600)
            addParameter(Parser,'EpochT',[0 150; 150 450; 450 750; 750 1050; 1050 1350])
            addParameter(Parser,'width',50)
            
            parse(Parser,obj,varargin{:})
            
            obj = Parser.Results.obj;
            dirs = Parser.Results.dirs;
            seqs = Parser.Results.seqs;
            EpochT = Parser.Results.EpochT;
            t = Parser.Results.t;
            width = Parser.Results.width;
            
            % Find mean eye speeds for each sequence
            for di = 1:length(dirs)
                R(:,:,:,di) = dynamicCohSeqConditionedRates(obj,...
                    'dirs',dirs(di),'width',width);
            end
            
            obj.Epoch.neurons = nan(max(diff(EpochT,1,2)),length(seqs),4,length(seqs)+1,length(dirs));
            for epochi = 1:size(EpochT,1)
                ttemp = 0:(EpochT(epochi,2)-EpochT(epochi,1));
                for controli = 1:(length(seqs)+1)
                    for seqi = 1:length(seqs)
                        for di = 1:length(dirs)
                            if controli-1 == 0
                                for uniti = 1:size(R,3)
                                    obj.Epoch.neurons(1:length(ttemp),seqi,epochi,controli,di,uniti) = ...
                                        R(t>=EpochT(epochi,1) & t<= EpochT(epochi,2),seqi,uniti,di);
                                    obj.Epoch.p(seqi,epochi,controli,di,uniti).neurons = ...
                                        pieceWiseExpRiseFitter(ttemp(:),...
                                            obj.Epoch.neurons(1:length(ttemp),seqi,epochi,controli,di,uniti),...
                                            [0,10,200,50],[Inf,Inf,Inf,200],[-Inf,-Inf,0,0]);
                                end
                            elseif controli-1 ~= seqi
                                obj.Epoch.neurons(1:length(ttemp),seqi,epochi,controli,di,:) = ...
                                    R(t>=EpochT(epochi,1) & t<= EpochT(epochi,2),seqi,:,di) -... 
                                    R(t>=EpochT(epochi,1) & t<= EpochT(epochi,2),controli-1,:,di);
                                
                                for uniti = 1:size(R,3)
                                    obj.Epoch.p(seqi,epochi,controli,di,uniti).neurons = ...
                                        pieceWiseExpRiseFitter(ttemp(:),...
                                            obj.Epoch.neurons(1:length(ttemp),seqi,epochi,controli,di,uniti),...
                                            [0,3,200,50],[Inf,Inf,Inf,200],[-Inf,-Inf,0,0]);
                                end
                            else
                                for uniti = 1:size(R,3)
                                    obj.Epoch.p(seqi,epochi,controli,di,uniti).neurons = nan(1,4);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function obj = findActive(obj,rateCutoff,cutWindow)
            if isempty(obj.R)
                obj.passCutoff = false(length(obj.unitIndex));
            else
                obj.passCutoff = (permute(max(obj.R(cutWindow(1):cutWindow(2),:,:,:),[],[1,2,4])*1000,[3,1,2]) - ...
                    permute(min(obj.R(cutWindow(1):cutWindow(2),:,:,:),[],[1,2,4])*1000,[3,1,2])) > rateCutoff;
            end
            obj.rateCutoff = rateCutoff;
            obj.cutWindow = cutWindow;
        end
        
        function obj = evaluatePreferredDirection(obj,win,varargin)
            % Parse inputs
            Parser = inputParser;
            
            addRequired(Parser,'obj')
            addRequired(Parser,'win')
            addParameter(Parser,'speeds',NaN)
            
            parse(Parser,obj,win,varargin{:})
            
            obj = Parser.Results.obj;
            win = Parser.Results.win;
            speeds = Parser.Results.speeds;
            
            % Find counts for each direction in window, marginalize speeds, coherences
            directions = unique(obj.directions);
            if isempty(directions)
                countsTotal = [];
            else
                for di = 1:length(directions)
                    counts = seqConditionedCounts(obj,'dirs',directions(di),'speeds',speeds,'win',win);
                    countsTotal(:,di) = permute(nansum(counts,[1,2]),[3,1,2]);
                end
            end
            [~,maxInds] = max(countsTotal,[],2);
            obj.preferredDirectionRelative = directions(maxInds);
            obj.preferredDirectionWin = win;
        end
        
        
        function C = findCovariance(obj,binT,speeds,dirs,...
                sequences,perturbations)
            res = nan(length(binT),1000,length(obj.unitIndex));
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for seqi = 1:length(sequences)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,NaN,sequences(seqi),perturbations);
                        for bi = 1:length(binT)
                            counts = obj.r(binT(bi) == obj.neuron_t,condLogical,:);
                            res(bi,ind:ind+sum(condLogical)-1,:) = counts*obj.filterWidth*2 - mean(counts*obj.filterWidth*2);
                        end
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(:,1:ind-1,:);
            for uniti = 1:length(obj.unitIndex)
                C(:,:,uniti) = cov(res(:,:,uniti)');
            end
            
        end
        
        function varCE = findVarCE(obj,speeds,dirs,sequences,perturbations)
            rtemp = obj.r*obj.filterWidth*2;
            res = nan(size(rtemp));
            res = permute(res,[1,3,2]);
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for seqi = 1:length(sequences)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,sequences(seqi),perturbations);
                        m(:,:,di,si,seqi) = mean(rtemp(:,condLogical,:),2);
                        res(:,:,ind:ind+sum(condLogical)-1) = permute(rtemp(:,condLogical,:),[1,3,2]) - ...
                            repmat(m(:,:,di,si,seqi),[1,1,sum(condLogical)]);
                        n(1,1,di,si,seqi) = sum(condLogical);
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(:,:,1:ind-1);
            M = sum(repmat(n,[size(m,1),size(m,2),1,1,1]).*m/sum(n(:)),[3,4,5]);
            V = var(res,[],3);
            FF = V./M;
            phi = min(FF,[],1);
            varCE = V - repmat(phi,[size(M,1),1]).*M;
        end
        
        
        function [rsc,pval] = spikeCountCorrelationWin(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'sequences',NaN)
            addParameter(Parser,'perturbations',NaN)
            addParameter(Parser,'win',[150,450])
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            dirs = Parser.Results.dirs;
            sequences = Parser.Results.sequences;
            perturbations = Parser.Results.perturbations;
            win = Parser.Results.win;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end            
            if isnan(sequences)
                squences = unique(obj.sequences);
            end               
            if isnan(perturbations)
                perturbations = unique(obj.perturbations);
            end
            
            % For each neuron and condition calculatate z-scores
            zscores = [];
            for seqi = 1:length(squences)
                for perti = 1:length(perturbations)
                    for di = 1:length(dirs)
                        
                        counts = conditionalCounts(obj,win,dirs(di),NaN,squences(seqi),...
                            perturbations(perti));
                        rstd = nanstd(counts,[],2);
                        rmean = nanmean(counts,2);
                        ztemp = (counts-repmat(rmean,[1,size(counts,2)]))./...
                            repmat(rstd,[1,size(counts,2)]);
                        zscores = cat(2,zscores,ztemp);
                    end
                end
            end
            
            % Find spike count correlations
            zscores = permute(zscores,[2,1]);
            [rsc, pval] = corrcoef(zscores);
        end
        
        
        function [cc, shufflecc] = conditionalCorrelograms(obj,unitsIndex,win,shuffleN,varargin)
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addRequired(Parser,'unitsIndex')
            addRequired(Parser,'win')
            addRequired(Parser,'shuffleN')
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'sequences',NaN)
            addParameter(Parser,'perturbations',NaN)
            
            parse(Parser,obj,unitsIndex,win,shuffleN,varargin{:})
            obj = Parser.Results.obj;
            win = Parser.Results.win;
            unitsIndex = Parser.Results.unitsIndex;
            shufffleN = Parser.Results.shuffleN;
            dirs = Parser.Results.dirs;
            sequences = Parser.Results.sequences;
            perturbations = Parser.Results.perturbations;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end            
            if isnan(sequences)
                speeds = unique(obj.sequences);
            end               
            if isnan(perturbations)
                perturbations = unique(obj.perturbations);
            end
            
            % Make edges from bin centers specified by win
            edges = [win-(win(2)-win(1))/2 win(end)+(win(2)-win(1))/2];
            
            cc = zeros([length(win),length(unitsIndex),length(unitsIndex)]);
            shufflecc = zeros([length(win),length(unitsIndex),length(unitsIndex),shuffleN]);
                
            condInds = trialSort(obj,dirs,NaN,NaN,NaN,sequences,perturbations);
            
            for triali = 1:length(condInds)
                spktimes = obj.spikeTimes{condInds(triali)}{1};
                clusters = obj.spikeTimes{condInds(triali)}{2};
                spktimes = spktimes(ismember(clusters,unitsIndex));
                clusters = clusters(ismember(clusters,unitsIndex));
                
                % Find differences and add to correlogram
                for uniti = 1:length(unitsIndex)
                    for unitj = 1:length(unitsIndex)
                        ts1 = spktimes(clusters == unitsIndex(uniti));
                        ts2 = spktimes(clusters == unitsIndex(unitj));
                        ind = 1;
                        for ti = 1:length(ts1)
                            d = ts1(ti) - ts2;
                            tdiffs = histc(d,edges);
                            cc(:,uniti,unitj) = cc(:,uniti,unitj) + tdiffs(1:end-1)';
                            
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
            h = figure;
            colors = colormap('lines');
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