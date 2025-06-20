function [dyn, gain, gainSTE] = dynamicCohBehavioralAnalysis(sname,varargin)
%%
%
%
%
%%

%% Defaults
dcp_default{1} = [];

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'sname')
addParameter(Parser,'dcp',dcp_default)
addParameter(Parser,'dcpObjectsFile',[])
addParameter(Parser,'trialList',1:5000)
addParameter(Parser,'speeds',NaN)
addParameter(Parser,'seqs',NaN)
addParameter(Parser,'win',[-Inf Inf])
addParameter(Parser,'pertTimes',150+[150:150:150*8])
addParameter(Parser,'pertMap',4:11)
addParameter(Parser,'motionChangeTimes',[150+0:300:1500])
addParameter(Parser,'forceRead',false)
addParameter(Parser,'dcpAccept',NaN)
addParameter(Parser,'outlierReject',false)
addParameter(Parser,'measurementMethod','peak2peak')
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'pertWin',300)

parse(Parser,sname,varargin{:})

sname = Parser.Results.sname;
dcp = Parser.Results.dcp;
dcpObjectsFile = Parser.Results.dcpObjectsFile;
trialList = Parser.Results.trialList;
speeds = Parser.Results.speeds;
seqs = Parser.Results.seqs;
win = Parser.Results.win;
pertTimes = Parser.Results.pertTimes;
pertMap = Parser.Results.pertMap;
motionChangeTimes = Parser.Results.motionChangeTimes;
forceRead = Parser.Results.forceRead;
dcpAccept = Parser.Results.dcpAccept;
outlierReject = Parser.Results.outlierReject;
measurementMethod = Parser.Results.measurementMethod;
directions = Parser.Results.directions;
pertWin = Parser.Results.pertWin;

%% Build dcp objects
if isempty(dcpObjectsFile) && isempty(dcp{1})
    % Build FileList
    potentialFiles = dir(['~/Projects/DynamicCoherencePhysiology/' sname]);
    regOut = regexpi({potentialFiles.name},'[0-9]{8}[a-z]{1,3}','match');
    ind = 0;
    for listi = 1:length(regOut)
        if ~isempty(regOut{listi}) && ~strcmp(regOut{listi}{1}(end-2:end),'plx') ...
                && ~strcmp(regOut{listi}{1}(end-1:end),'kk') ...
                && ~any(strcmp(regOut{listi}{1},excludeList))
            ind = ind+1;
            FileList{ind} = regOut{listi}{1};
        end
    end
    dcp = dcpPrelim(subject,FileList,extractSpikes);
elseif isempty(dcp{1})
    load(['~/Projects/DynamicCoherencePhysiology/' sname '/dcpObjects/' dcpObjectsFile],'dcp')
    if ~isnan(dcpAccept)
        dcp = dcp(dcpAccept);
    end
end

%% Get data from dynamicCoh experiments
ind = 0;
emptyT = true;
while emptyT
    ind = ind+1;
    dynCoh = dynamicCohObj(dcp{ind}.sname,dcp{ind}.datapath);
    dynCoh = dynamicCohTrials(dynCoh,trialList,forceRead);
    dynCoh = addCoh(dynCoh);
    
    dyn.t = dynCoh.eye_t';
    dyn.coh = dynCoh.coh;
    emptyT = isempty(dyn.t);
end
dyn.conditions.directions = nan(max(trialList)*length(dcp),1);
dyn.conditions.speeds = nan(max(trialList)*length(dcp),1);
dyn.conditions.seq = nan(max(trialList)*length(dcp),1);
dyn.conditions.perts = nan(max(trialList)*length(dcp),1);
dyn.eye.hvel = nan(length(dyn.t),max(trialList)*length(dcp));
dyn.eye.vvel = nan(length(dyn.t),max(trialList)*length(dcp));
ind = 1;
for filei = 1:length(dcp)
    
    dynCoh = dynamicCohObj(dcp{filei}.sname,dcp{filei}.datapath);
    dynCoh = dynamicCohTrials(dynCoh,trialList,forceRead);
    
    if ~isempty(dynCoh.trialtype)
        dyn.conditions.directions(ind:ind+length(dynCoh.directions)-1) = dynCoh.directions;
        dyn.conditions.speeds(ind:ind+length(dynCoh.directions)-1) = dynCoh.speeds;
        dyn.conditions.seq(ind:ind+length(dynCoh.directions)-1) = dynCoh.sequences;
        dyn.conditions.perts(ind:ind+length(dynCoh.directions)-1) = dynCoh.perturbations;
        
        dyn.eye.hvel(:,ind:ind+length(dynCoh.directions)-1) = vertcat(dynCoh.eye(:).hvel)';
        dyn.eye.vvel(:,ind:ind+length(dynCoh.directions)-1) = vertcat(dynCoh.eye(:).vvel)';
        
        ind = ind+length(dynCoh.directions);
    end
end
dyn.conditions.directions = dyn.conditions.directions(1:ind-1);
dyn.conditions.speeds = dyn.conditions.speeds(1:ind-1);
dyn.conditions.seq = dyn.conditions.seq(1:ind-1);
dyn.conditions.perts = dyn.conditions.perts(1:ind-1);
dyn.eye.hvel = dyn.eye.hvel(:,1:ind-1);
dyn.eye.vvel = dyn.eye.vvel(:,1:ind-1);
dyn.eye.speed = sqrt(dyn.eye.hvel.^2 + dyn.eye.vvel.^2);

%% Find conditional means, marginalize direction and perturbations
if isnan(speeds)
    speeds = unique(dyn.conditions.speeds);
end
if isnan(seqs)
    seqs = unique(dyn.conditions.seq);
end
for si = 1:length(speeds)
    for seqi = 1:length(seqs)
        acceptvec = dyn.conditions.speeds == speeds(si) & dyn.conditions.seq == seqs(seqi) & ...
            ismember(dyn.conditions.directions,directions);
        dyn.eye.mean(:,si,seqi) = nanmean(dyn.eye.speed(:,acceptvec),2);
        dyn.eye.ste(:,si,seqi) = nanstd(dyn.eye.speed(:,acceptvec),[],2)./sqrt(sum(acceptvec));
    end
end

%% Find perturbation response, marginalize direction and speed (normalized)
perts = unique(dyn.conditions.perts);
for seqi = 1:length(seqs)
    stemp = [];
        
        for pi = 1:length(perts)
            acceptvec = dyn.conditions.seq == seqs(seqi) & ...
                dyn.conditions.perts == perts(pi) &...
                ismember(dyn.conditions.directions,directions);
            dyn.eye.pert.m(:,seqi,pi) = nanmean(dyn.eye.speed( dyn.t>=win(1) & dyn.t<=win(2),acceptvec)...
                ./repmat(dyn.conditions.speeds(acceptvec)',[length(dyn.t>=win(1) & dyn.t<=win(2)),1]),2);
            dyn.eye.pert.ste(:,seqi,pi) = nanstd(dyn.eye.speed( dyn.t>=win(1) & dyn.t<=win(2),acceptvec)...
                ./repmat(dyn.conditions.speeds(acceptvec)',[length(dyn.t>=win(1) & dyn.t<=win(2)),1]),[],2)/sqrt(sum(acceptvec));
            dyn.eye.pert.N(:,seqi,pi) = sum(acceptvec);
            
            if perts(pi) == 0
                dyn.eye.pert.t(seqi,pi) = NaN;
                dyn.eye.pert.res(seqi,pi) = 0;
                dyn.eye.pert.coh(seqi,pi) = NaN;
                dyn.eye.pert.resSTE(seqi,pi) = 0;
                dyn.eye.pert.resControl(seqi,pi) = 0;
                dyn.eye.pert.resControlSTE(seqi,pi) = 0;
            else
                dyn.eye.pert.t(seqi,pi) = pertTimes(perts(pi) == pertMap);
                if any(perts == 0)
                    diffTemp = dyn.eye.pert.m(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin,seqi,pi) - ...
                        dyn.eye.pert.m(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin,seqi,perts == 0);
                    steTemp = sqrt(dyn.eye.pert.ste(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin,seqi,pi).^2 + ...
                        dyn.eye.pert.ste(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin,seqi,perts == 0).^2);
                    
                    [maxRes, maxInd] = max(diffTemp);
                    [minRes, minInd] = min(diffTemp);
                    dyn.eye.pert.res(seqi,pi) = maxRes - minRes;
                    dyn.eye.pert.resSTE(seqi,pi) = sqrt(steTemp(maxInd).^2 + steTemp(minInd).^2);
                    dyn.eye.pert.coh(seqi,pi) = dyn.coh(dyn.t == dyn.eye.pert.t(seqi,pi),seqi);
                    
                    diffTemp2 = dyn.eye.pert.m(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi)-pertWin & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi),seqi,pi) - ...
                        dyn.eye.pert.m(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi)-pertWin & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi),seqi,perts == 0);                    
                    steTemp2 = sqrt(dyn.eye.pert.ste(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi)-pertWin & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi),seqi,pi).^2 + ...
                        dyn.eye.pert.ste(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi)-pertWin & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi),seqi,perts == 0).^2);
                    
                    [maxRes, maxInd] = max(diffTemp2);
                    [minRes, minInd] = min(diffTemp2);
                    dyn.eye.pert.resControl(seqi,pi) = maxRes - minRes;
                    dyn.eye.pert.resControlSTE(seqi,pi) = sqrt(steTemp2(maxInd).^2 + steTemp2(minInd).^2);
                    
                else
                    acceptvec = dyn.conditions.seq == seqs(seqi) & ...
                        dyn.conditions.perts ~= perts(pi);
                    control = nanmean(dyn.eye.speed( dyn.t>=win(1) & dyn.t<=win(2),acceptvec)...
                        ./repmat(dyn.conditions.speeds(acceptvec)',[length(dyn.t>=win(1) & dyn.t<=win(2)),1]),2);
                    controlSTE = nanstd(dyn.eye.speed( dyn.t>=win(1) & dyn.t<=win(2),acceptvec)...
                        ./repmat(dyn.conditions.speeds(acceptvec)',[length(dyn.t>=win(1) & dyn.t<=win(2)),1]),2)./sqrt(sum(acceptvec));
                    diffTemp = dyn.eye.pert.m(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin,seqi,pi) - ...
                        control(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin);                    
                    steTemp = sqrt(dyn.eye.pert.ste(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin,seqi,pi).^2 + ...
                        controlSTE(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)+pertWin).^2);
                    
                    [maxRes, maxInd] = max(diffTemp);
                    [minRes, minInd] = min(diffTemp);
                    dyn.eye.pert.res(seqi,pi) = maxRes - minRes;
                    dyn.eye.pert.resSTE(seqi,pi) = sqrt(steTemp(maxInd).^2 + steTemp(minInd).^2);
                    dyn.eye.pert.coh(seqi,pi) = dyn.coh(dyn.t == dyn.eye.pert.t(seqi,pi),seqi);
                    
                    diffTemp2 = dyn.eye.pert.m(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)-pertWin,seqi,pi) - ...
                        control(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)-pertWin);                    
                    steTemp2 = sqrt(dyn.eye.pert.ste(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)-pertWin,seqi,pi).^2 + ...
                        control(...
                        dyn.t >= dyn.eye.pert.t(seqi,pi) & ...
                        dyn.t <= dyn.eye.pert.t(seqi,pi)-pertWin).^2);
                    
                    [maxRes, maxInd] = max(diffTemp2);
                    [minRes, minInd] = min(diffTemp2);
                    dyn.eye.pert.resControl(seqi,pi) = maxRes - minRes;
                    dyn.eye.pert.resControlSTE(seqi,pi) = sqrt(steTemp2(maxInd).^2 + steTemp2(minInd).^2);
                end
            end
            
            if outlierReject 
                for iter = 1:2
                    mtemp = nanmean(dyn.eye.init{si,seqi});
                    stdtemp = nanstd(dyn.eye.init{si,seqi});
                    dyn.eye.init{si,seqi} = dyn.eye.init{si,seqi}(...
                        dyn.eye.init{si,seqi} < mtemp + 1.96*stdtemp & ...
                        dyn.eye.init{si,seqi} > mtemp - 1.96*stdtemp);
                end
            end
        end
    
end
% 
% [gain.B(:,seqi), gain.BINT(:,seqi), ~, ~, gain.STATS(seqi)] = ...
%         regress(horzcat(dyn.eye.init{:,seqi})',[stemp ones(size(stemp))]);

%% Measure gain
for seqi = 1:length(seqs)
    for perti = 2:length(perts)
        switch measurementMethod
            case 'integral'
                temp = (dyn.eye.pert.m(dyn.eye.pert.t(seqi,perti):dyn.eye.pert.t(seqi,perti)+300,seqi,perti) - ...
                    dyn.eye.pert.m(dyn.eye.pert.t(seqi,perti):dyn.eye.pert.t(seqi,perti)+300,seqi,1));
                gain(seqi,perti) = sum( abs(temp-mean(temp) ))/8;                   % 8 is the integral of the absolute value of a full cycle perturbation with amplitude 2 deg/s
        
            case 'peak2peak'
                gain(seqi,perti) = (dyn.eye.pert.res(seqi,perti)-dyn.eye.pert.resControl(seqi,perti))/0.4;                         % peak to peak amplitude of perturbation was 40% of target speed)
                gainSTE(seqi,perti) = sqrt(dyn.eye.pert.resSTE(seqi,perti).^2 + dyn.eye.pert.resControlSTE(seqi,perti).^2)/0.4;
        end
    end
end

%% Plotting

%% Mean dynamicCoh response
figure;
colors = colormap('lines');
close(gcf)
figure('Name','Mean dynCoh response')
for seqi = 1:length(seqs)
    for si = 1:length(speeds)
        subplot(1,length(speeds),si)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = colors(seqi,:);
        myPatch(dyn.t(:),dyn.eye.mean(:,si,seqi),dyn.eye.ste(:,si,seqi),'patchProperties',patchProps);
        hold on
        
        plot(dyn.t,dyn.eye.mean(:,si,seqi),'Color',colors(seqi,:),'LineWidth',2)
        
        xlabel('Time from motion onset (ms)')
        ylabel('Eye speeds (deg/s)')
        title(['Target speed = ' num2str(speeds(si)) ' (deg/s)']);
        plotVertical(motionChangeTimes);
    end
end
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

%% Change in eye speeds
figure;
colors = colormap('lines');
close(gcf)
figure('Name','Mean dynCoh response')
subplot(2,1,1)
for seqi = 1:length(seqs)
    for si = 1:length(speeds)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = colors(seqi,:);
        myPatch(dyn.t(:),dyn.eye.mean(:,si,seqi)-dyn.eye.mean(:,si,5),sqrt(dyn.eye.ste(:,si,seqi).^2+dyn.eye.ste(:,si,5).^2),'patchProperties',patchProps);
        hold on
        
        plot(dyn.t,dyn.eye.mean(:,si,seqi)-dyn.eye.mean(:,si,5),'Color',colors(seqi,:),'LineWidth',2)
        
        xlabel('Time from motion onset (ms)')
        ylabel('Eye speeds (deg/s)')
        title(['Target speed = ' num2str(speeds(si)) ' (deg/s)']);
        plotVertical(motionChangeTimes);
    end
end
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

% Collate across identical conditions up to 750 ms
subplot(2,1,2)
seqsToCombine = [1,3; 2,4];
for seqi = 1:size(seqsToCombine,1)
    for si = 1:length(speeds)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = colors(seqi,:);
        
        acceptvec = dyn.conditions.speeds == speeds(si) & ismember(dyn.conditions.seq,seqsToCombine(seqi,:)) & ...
            ismember(dyn.conditions.directions,directions);
        meyespeeds = nanmean(dyn.eye.speed(:,acceptvec),2);
        steeyespeeds = std(dyn.eye.speed(:,acceptvec),[],2)/sqrt(sum(acceptvec));
                
        myPatch(dyn.t(:),meyespeeds-dyn.eye.mean(:,si,5),sqrt(steeyespeeds.^2+dyn.eye.ste(:,si,5).^2),'patchProperties',patchProps);
        hold on
        
        plot(dyn.t,meyespeeds-dyn.eye.mean(:,si,5),'Color',colors(seqi,:),'LineWidth',2)
        
        xlabel('Time from motion onset (ms)')
        ylabel('Eye speeds (deg/s)')
        title(['Target speed = ' num2str(speeds(si)) ' (deg/s)']);
        plotVertical(motionChangeTimes);
    end
end
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])


%% Perturbations
figure('Name','Pertubation response','Position',[855 54 583 1242])
pertInds = 1:length(perts);
pertInds = pertInds(perts ~= 0);
ind = 0;
for pi = pertInds
    ind = ind+1;
    subplot(sum(perts~=0),1,ind)
    for seqi = 1:length(seqs)
        if any(perts == 0)
            plot(dyn.t,dyn.eye.pert.m(:,seqi,perts == perts(pi)) - dyn.eye.pert.m(:,seqi,perts == 0),'Color',colors(seqi,:))
        else
            acceptvec = dyn.conditions.seq == seqs(seqi) & ...
                dyn.conditions.perts ~= perts(pi);
            control = nanmean(dyn.eye.speed( dyn.t>=win(1) & dyn.t<=win(2),acceptvec)...
                ./repmat(dyn.conditions.speeds(acceptvec)',[length(dyn.t>=win(1) & dyn.t<=win(2)),1]),2);
            
            plot(dyn.t,dyn.eye.pert.m(:,seqi,perts == perts(pi)) - control,'Color',colors(seqi,:))
        end
        hold on
        
        if perts(pi) == 0
            legendEntry{seqi} = 'Control';
        else
            legendEntry{seqi} = num2str(dyn.eye.pert.coh(seqi,pi));
        end
    end
    ylim(0.2*[-1 1])
    plotVertical(motionChangeTimes);
    if perts(pi) ~= 0
        lineProps.Color = [1 0 0];
        plotVertical(pertTimes(pertMap == perts(pi)),'lineProperties',lineProps);
    end
    legend(legendEntry)
    xlabel('Time from motion onset (ms)')
    ylabel('Perturbation response / target speed')
end
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

figure('Name','Feedforward gain estimate')
symbols = {'d','o','s'};
for pi = 1:3
    for seqi = 1:length(seqs)
        %     plot(dyn.eye.pert.coh(seqi,:)',dyn.eye.pert.res(seqi,:)'/0.2/2,...
        %         'o-','Color',colors(seqi,:))
        errorbar(dyn.eye.pert.coh(seqi,pi+1)',gain(seqi,pi+1)',gainSTE(seqi,pi+1)*1.64,[],...
            symbols{pi},'Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),'MarkerSize',8)
        hold on
    end
end
plotHorizontal(0);
xlabel('Coherence')
ylabel('Gain')
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

figure('Name','Feedforward gain estimate')
symbols = {'d','o','s'};
for pi = 1:3
    for seqi = 1:length(seqs)
        %     plot(dyn.eye.pert.coh(seqi,:)',dyn.eye.pert.res(seqi,:)'/0.2/2,...
        %         'o-','Color',colors(seqi,:))
        errorbar(dyn.eye.pert.coh(seqi,pi+1)',gain(seqi,pi+1)',gainSTE(seqi,pi+1)*1.64,[],...
            symbols{pi},'Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),'MarkerSize',8)
        hold on
    end
end
plotHorizontal(0);
xlabel('Coherence')
ylabel('Gain')
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

