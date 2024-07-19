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
        objType = 'initiateCohObj';
        coh;
        cohTuning;
        r;
        R;
        Rste;
        
        preferredDirection;
        preferredDirectionRelative;
        
        saveLocation;
    end
    
    properties (SetAccess = private)
        filterWidth;
        preferredDirectionWin;
        rateCutoff;
        cutWindow;
        passCutoff;
        
        reward;
        rotation;
        
        eye_t;
        neuron_t;
    end
    
    methods
        %% Core methods
        function obj = initiateCohObj(sname,datapath)
        %Constructor
            obj = obj@dcpObj(sname,datapath);
            obj.cohTuning.parameters = nan(1,3);
        end 
        
        %% initiateCohTrials
        function obj = initiateCohTrials(obj,trialNs,acceptInitCohPert)
        % Add direction preference trials to data object
            if ~exist('acceptInitCohPert','var')
                acceptInitCohPert = false;
            end
        
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
                        %disp([obj.datapath '/' files(ti+fileInx-1).name])
                        file = readcxdata([obj.datapath '/' files(ti+fileInx-1).name]);
                        
                        trialname = file.trialname;
                        if ~isempty(trialname) && acceptInitCohPert && strcmp(trialname(1:7),'initCoh')
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
                            
                            obj.rotation(ind,1) = file.key.iVelTheta/1000;
                            
                            obj.reward(ind,1) = file.key.iRewLen1;
                            
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
                            
                        elseif ~isempty(trialname) && strcmp(trialname(1:7),'initCoh') && ~strcmp(trialname(1:11),'initCohPert')
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
                            
                            obj.rotation(ind,1) = file.key.iVelTheta/1000;
                            
                            obj.reward(ind,1) = file.key.iRewLen1;
                            
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
        
        %% Behavioral analysis methods
        
        function [C, res] = findCovarianceBehavior(obj,speeds,cohs,dirs,win,interpolationMethod,interpolationThreshold,detectPoorPursuitThreshold)
            if ~exist('interpolationMethod','var')
                interpolationMethod = 'Linear';
            end
            if ~exist('interpolationThreshold','var')
                interpolationThreshold = 0;
            end
            if ~exist('detectPoorPursuitThreshold','var')
                detectPoorPursuitThreshold = Inf;
            end
            res = nan(1000,length(obj.eye_t(obj.eye_t>=win(1) & obj.eye_t<= win(2))));
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        eh = vertcat(obj.eye(:).hvel);
                        ev = vertcat(obj.eye(:).vvel);
                        eSpeed = sqrt(eh(condLogical,:).^2 + ...
                            ev(condLogical,:).^2);
                        
                        rmse = sqrt(nanmean(...
                            (eSpeed(:,obj.eye_t > 200 & obj.eye_t < 1200)-speeds(si)).^2,2))...
                            ./speeds(si);
                        
                        meetsThreshold = sum(~isnan(eSpeed),2) > interpolationThreshold & rmse < detectPoorPursuitThreshold;
                        eSpeed = eSpeed(meetsThreshold,:);
                        
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
                        
                        res(ind:ind+size(eSpeed,1)-1,:) = ...
                            eSpeed(:,obj.eye_t>=win(1) & obj.eye_t<= win(2)) - ...
                            mean(eSpeed(:,obj.eye_t>=win(1) & obj.eye_t<= win(2)),1);
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(1:ind-1,:);
            stdRes = std(res(:));
            res = res(~any(abs(res) > 1.5*stdRes,2),:);
            C = cov(res);
            
        end
        
        %% Neural analysis methods
        function [r,rste] = conditionalRates(obj,width,directions,speeds,...
                cohs,t)
            if ~exist('t','var')
                t = -100:1600;
            end
            rAll = obj.calcRates(width,'t',t);
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,cohs);
            r = mean(rAll(:,condLogical,:),2);
            rste = sqrt(var(rAll(:,condLogical,:),[],2)/sum(condLogical));
        end
        
        function obj = cohConditionedRates(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'width',50)
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'cohs',NaN)
            addParameter(Parser,'t',-100:1600)
            addParameter(Parser,'marginalizeDirection',true)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            width = Parser.Results.width;
            dirs = Parser.Results.dirs;
            speeds = Parser.Results.speeds;
            locs = Parser.Results.locs;
            cohs = Parser.Results.cohs;
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
            if isnan(cohs)
                cohs = unique(obj.coh);
            end
            
            % Condition on sequences
            for speedi = 1:length(speeds)
                for cohi = 1:length(cohs)
                    if marginalizeDirection
                        [R(:,speedi,cohi,:),Rste(:,speedi,cohi,:)] = ...
                            conditionalRates(obj,width,...
                            dirs,speeds(speedi),cohs(cohi),t);
                    else
                        for di = 1:length(dirs)
                            [R(:,speedi,cohi,:,di),Rste(:,speedi,cohi,:,di)] = ...
                                conditionalRates(obj,width,...
                                dirs(di),speeds(speedi),cohs(cohi),t);
                        end
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
                cohs)
            countsAll = obj.countSpikes(win);
            if exist('countsAll','var')
                [~,condLogical] = trialSort(obj,directions,speeds,NaN,cohs);
                counts = countsAll(:,condLogical);
            else
                [~,condLogical] = trialSort(obj,directions,speeds,NaN,cohs);
                counts = zeros(1,sum(condLogical));
            end
        end
        
        function counts = cohConditionedCounts(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'win',[176,225])
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'cohs',NaN)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            win = Parser.Results.win;
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
            counts = nan([1000,length(speeds),length(cohs),length(obj.unitIndex)]);
            maxtrials = 0;
            for speedi = 1:length(speeds)
                for cohi = 1:length(cohs)
                    countsTemp = permute(conditionalCounts(obj,win,...
                        dirs,speeds(speedi),cohs(cohi)),[4,2,3,1]);
                    if size(countsTemp,1) > 1
                        counts(1:size(countsTemp,1),speedi,cohi,:) = countsTemp;
                        maxtrials = max([maxtrials size(countsTemp,1)]);
                    else
                        counts(1:length(countsTemp),speedi,cohi,:) = countsTemp;
                        maxtrials = max([maxtrials length(countsTemp)]);
                    end
                end
            end
            counts = counts(1:maxtrials,:,:,:);
        end
               
        function obj = findActive(obj,rateCutoff,cutWindow)
            if isempty(obj.R)
                obj.passCutoff = false(length(obj.unitIndex));
            else
                obj.passCutoff = (permute(max(obj.R(cutWindow(1):cutWindow(2),:,:,:,:),[],[1,2,3,5])*1000,[4,1,2,3]) - ...
                    permute(min(obj.R(cutWindow(1):cutWindow(2),:,:,:,:),[],[1,2,3,5])*1000,[4,1,2,3])) > rateCutoff;
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
                    counts = cohConditionedCounts(obj,'dirs',directions(di),'speeds',speeds,'win',win);
                    countsTotal(:,di) = permute(nansum(counts,[1,2,3]),[4,1,2,3]);
                end
            end
            [~,maxInds] = max(countsTotal,[],2);
            obj.preferredDirectionRelative = directions(maxInds);
            obj.preferredDirectionWin = win;
        end
        
        function C = findCovariance(obj,binT,speeds,cohs,dirs)
            res = nan(length(binT),1000,length(obj.unitIndex));
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        for bi = 1:length(binT)
                            counts = obj.r(binT(bi) == obj.neuron_t,condLogical,:);
                            res(bi,ind:ind+sum(condLogical)-1,:) = counts*obj.filterWidth*2 - mean(counts*obj.filterWidth*2);
%                             if length(obj.unitIndex) == 1
%                                 resTemp(bi,:) = counts*obj.filterWidth*2 - mean(counts*obj.filterWidth*2);
%                             else
%                                 resTemp(bi,:,:) = counts*obj.filterWidth*2 - mean(counts*obj.filterWidth*2);
%                             end
                        end
%                         res(:,ind:ind+size(resTemp,2)-1,:) = resTemp;
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(:,1:ind-1,:);
            for uniti = 1:length(obj.unitIndex)
                C(:,:,uniti) = cov(res(:,:,uniti)');
            end
            
        end
        
        function varCE = findVarCE(obj,speeds,cohs,dirs)
            rtemp = obj.r*obj.filterWidth*2;
            res = nan(size(rtemp));
            res = permute(res,[1,3,2]);
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        m(:,:,di,si,ci) = mean(rtemp(:,condLogical,:),2);
                        res(:,:,ind:ind+sum(condLogical)-1) = permute(rtemp(:,condLogical,:),[1,3,2]) - ...
                            repmat(m(:,:,di,si,ci),[1,1,sum(condLogical)]);
                        n(1,1,di,si,ci) = sum(condLogical);
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
        
        function [res, m, n] = calcResiduals(obj,speeds,cohs,dirs)
            rtemp = obj.r*obj.filterWidth*2;
            res = nan(size(rtemp));
            res = permute(res,[1,3,2]);
            ind = 1;
            for di = 1:length(dirs)
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        m(:,:,di,si,ci) = mean(rtemp(:,condLogical,:),2);
                        res(:,:,ind:ind+sum(condLogical)-1) = permute(rtemp(:,condLogical,:),[1,3,2]) - ...
                            repmat(m(:,:,di,si,ci),[1,1,sum(condLogical)]);
                        n(1,1,di,si,ci) = sum(condLogical);
                        ind = ind+sum(condLogical);
                    end
                end
            end
            res = res(:,:,1:ind-1);
        end
        
        function [rsc,pval] = spikeCountCorrelation(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'cohs',NaN)
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
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
            if isnan(cohs)
                cohs = unique(obj.coh);
            end
            
            % For each neuron and condition calculatate z-scores
            zscores = [];
            for si = 1:length(speeds)
                for ci = 1:length(cohs)
                    for di = 1:length(dirs)
                        [~,condLogical] = trialSort(obj,dirs(di),speeds(si),NaN,cohs(ci));
                        r = obj.r(:,condLogical,:);
                        rstd = nanstd(r,[],2);
                        rmean = nanmean(r,2);
                        ztemp = (r-repmat(rmean,[1,sum(condLogical),1]))./...
                            repmat(rstd,[1,sum(condLogical),1]);
                        zscores = cat(2,zscores,ztemp);
                    end
                end
            end
            
            % Find spike count correlations
            zscores = permute(zscores,[2,3,1]);
            for tind = 1:size(zscores,3)
                [rsc(:,:,tind),pval(:,:,tind)] = corrcoef(zscores(:,:,tind));
            end
        end
        
        
        function [rsc,pval] = spikeCountCorrelationWin(obj,varargin)
            
            % Parse inputs
            Parser = inputParser;
            addRequired(Parser,'obj')
            addParameter(Parser,'dirs',NaN)
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'cohs',NaN)
            addParameter(Parser,'win',[150,450])
            
            parse(Parser,obj,varargin{:})
            obj = Parser.Results.obj;
            dirs = Parser.Results.dirs;
            speeds = Parser.Results.speeds;
            locs = Parser.Results.locs;
            cohs = Parser.Results.cohs;
            win = Parser.Results.win;
            
            if isnan(dirs)
                dirs = unique(obj.directions);
            end            
            if isnan(speeds)
                speeds = unique(obj.speeds);
            end               
            if isnan(cohs)
                cohs = unique(obj.coh);
            end
            
            % For each neuron and condition calculatate z-scores
            zscores = [];
            for si = 1:length(speeds)
                for ci = 1:length(cohs)
                    for di = 1:length(dirs)
                        
                        counts = conditionalCounts(obj,win,dirs(di),speeds(si),...
                            cohs(ci));
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
            addParameter(Parser,'speeds',NaN)
            addParameter(Parser,'locs',NaN)
            addParameter(Parser,'cohs',NaN)
            
            parse(Parser,obj,unitsIndex,win,shuffleN,varargin{:})
            obj = Parser.Results.obj;
            win = Parser.Results.win;
            unitsIndex = Parser.Results.unitsIndex;
            shufffleN = Parser.Results.shuffleN;
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
            if isnan(cohs)
                cohs = unique(obj.coh);
            end
            
            % Make edges from bin centers specified by win
            edges = [win-(win(2)-win(1))/2 win(end)+(win(2)-win(1))/2];
            
            cc = zeros([length(win),length(unitsIndex),length(unitsIndex)]);
            shufflecc = zeros([length(win),length(unitsIndex),length(unitsIndex),shuffleN]);
                
            condInds = trialSort(obj,dirs,speeds,NaN,cohs);
            
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
        
        function [W, X, deltaT] = ...
                eyeRegressionModel(obj,speed,dir,coh,win,smoothWin,delays)
            
            % Find average eye position, velocity, and acceleration for
            % condition
            acceptVec = obj.speeds == speed & obj.directions == dir & obj.coh == coh;
            positions = vertcat(obj.eye(acceptVec).hpos);
            position = nanmean(positions,1);
            velocities = vertcat(obj.eye(acceptVec).hvel);
            velocity = nanmean(velocities,1);
            acceleration = [0 diff(velocity)];  % TODO
            
            velocity = smooth(velocity,smoothWin);
            acceleration = smooth(acceleration,smoothWin);
            
            X = [position' velocity acceleration ones(size(position'))];
            
            tVec = obj.eye_t >= win(1) & obj.eye_t <= win(2);
            X = X(tVec,:);
            
            % Regress with average firing rate for that condition
            speedInd = find(speed == unique(obj.speeds));
            cohInd = find(coh == unique(obj.coh));
            dirInd = find(dir == unique(obj.directions));
            for neuroni = 1:size(obj.R,4)
                FR = obj.R(:,speedInd,cohInd,neuroni,dirInd);
                
                for delayi = 1:length(delays)
                    tVec2 = obj.neuron_t >= win(1)+delays(delayi) & ...
                        obj.neuron_t <= win(2)+delays(delayi);
                    
                    [Wall(:,neuroni,delayi),~,res] = regress(FR(tVec2),X);
                    sse(neuroni,delayi) = sum(res.^2);
                end
                
                [~,minInd] = min(sse(neuroni,:),[],2);
                W(:,neuroni) = Wall(:,neuroni,minInd);
                deltaT(neuroni) = delays(minInd);
            end
        end
                
        %% Plotting methods
        function [h,sh,colors] = initiateCohMeanEye(obj,dirs,normalize,window,figureType)
        % Plots mean eye speed for each sequence and target speed
            if ~exist('normalize','var')
                normalize = true;
            end
            if ~exist('window','var')
                window = [0 300];
            end
            if ~exist('figureType','var')
                figureType = 'matlab';
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
                    switch figureType
                        case 'fyp'
                            plotMeanEyeSpeed(obj,condLogical,'normalizer',nrm,...
                                'h',h,'sh',sh(hi),'color',colors(si,:),'plotPatch',false);
                        case 'matlab'
                            plotMeanEyeSpeed(obj,condLogical,'normalizer',nrm,...
                                'h',h,'sh',sh(hi),'color',colors(si,:),'plotPatch',true);
                        otherwise
                            error(['Figure type ' figureType ' not recognized!'])
                    end
                end
                axis tight
                ax(hi,:) = axis;
            end
            for hi = 1:length(cohs)
                subplot(sh(hi))
%                 axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
                axis([window min(ax(:,3)) max(ax(:,4))])
                plotVertical(150);
                xlabel('Time from motion onset (ms)')
                if normalize
                    plotHorizontal(1);
                    switch figureType
                        case 'fyp'
                            ylabel('Eye speed / Target speed')
                            title([num2str(cohs(hi)) '% coherence'])
                        case 'matlab'
                            ylabel('$\frac{\textrm{Eye speed}}{\textrm{Target speed}}$')
                            mymakeaxis(gca,'xytitle',[num2str(cohs(hi)) '\% \textrm{coherence}'],...
                                'interpreter','latex','yticks',[0.1 1])
                        otherwise
                            error(['Figure type ' figureType ' not recognized!'])
                    end
                else
                    switch figureType
                        case 'fyp'
                            ylabel('Eye speed (deg/s)')
                            title([num2str(cohs(hi)) '% coherence'])
                        case 'matlab'
                            ylabel('Eye speed (deg/s)')
                            mymakeaxis(gca,'xytitle',[num2str(cohs(hi)) '\% \textrm{coherence}'],...
                                'interpreter','latex')
                        otherwise
                            error(['Figure type ' figureType ' not recognized!'])
                    end                  
                end
            end
            
            switch figureType
                case 'fyp'
                    for si = 1:length(speeds)
                        leg{si} = [num2str(speeds(si)) ' deg/s'];
                    end
                case 'matlab'
                    for li = 1:(2*length(speeds)+2)
                        leg{li} = '';
                    end
                    for si = 1:length(speeds)
                        leg{2*si} = [num2str(speeds(si)) ' deg/s'];
                    end
            end
            legend(leg)
        end
    end
end