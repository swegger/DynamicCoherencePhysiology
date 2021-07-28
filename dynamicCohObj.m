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
        rateCutoff;
        passCutoff;
        
        Epoch;
        
        saveLocation;
    end
    
    properties (SetAccess = private)
        filterWidth;
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
        
        %% addCoh
        function obj = addCoh(obj)
            t = -100:1600;
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
        
        %% Neural analysis methods
        function [r,rste] = conditionalRates(obj,width,directions,speeds,...
                sequences,perturbations,t)
            if ~exist('t','var')
                t = -100:1600;
            end
            
            % Determine if any spikies were recorded
            units = [];
            for triali = 1:length(obj.spikeTimes)
                units = [units obj.spikeTimes{triali}{2}];
            end
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                sequences,perturbations);
            if isempty(units)
                warning(['No spike times recorded for dynamicCohObj associated with ' obj.datapath(end-8:end)])
                r = nan(length(t),sum(condLogical));
                rste = nan(length(t),sum(condLogical));
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
                        dirs,speeds,seqs(seqi),perts);
                else
                    for di = 1:length(dirs)
                        [R(:,seqi,:,di),Rste(:,seqi,:,di)] = conditionalRates(obj,width,...
                            dirs(di),speeds,seqs(seqi),perts);                    
                    end
                end
            end
            obj.R = R;
            obj.Rste = Rste;
            obj.filterWidth = width;
            
        end
        
        function counts = conditionalCounts(obj,win,directions,speeds,...
                sequences,perturbations)
            countsAll = obj.countSpikes(win);
            [~,condLogical] = trialSort(obj,directions,speeds,NaN,NaN,...
                sequences,perturbations);
            counts = countsAll(:,condLogical);
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
        
        function obj = findActive(obj)
            obj.passCutoff = permute(max(obj.R,[],[1,2,4])*1000,[3,1,2])>obj.rateCutoff;
        end
        
        function obj = set.rateCutoff(obj,rateCutoff)
            obj.rateCutoff = rateCutoff;
            obj = findActive(obj);
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