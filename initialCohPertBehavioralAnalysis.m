function [init, gain] = initialCohPertBehavioralAnalysis(sname,varargin)
%%
%
%
%
%%

%% Defaults
plotSamps_default.On = true;
plotSamps_default.n = 100;
plotSamps_default.replacement = false;
dcp_default{1} = [];

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'sname')
addParameter(Parser,'dcp',dcp_default)
addParameter(Parser,'dcpObjectsFile',[])
addParameter(Parser,'trialList',1:5000)
addParameter(Parser,'win',[150 200])
addParameter(Parser,'plotSamps',plotSamps_default)
addParameter(Parser,'outlierReject',false)
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'keep_pert3always0deg',false)
addParameter(Parser,'pertWin',300)

parse(Parser,sname,varargin{:})

sname = Parser.Results.sname;
dcp = Parser.Results.dcp;
dcpObjectsFile = Parser.Results.dcpObjectsFile;
trialList = Parser.Results.trialList;
win = Parser.Results.win;
plotSamps = Parser.Results.plotSamps;
outlierReject = Parser.Results.outlierReject;
directions = Parser.Results.directions;
keep_pert3always0deg = Parser.Results.keep_pert3always0deg;
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
end

pert3always0deg = false(1,length(dcp));
for filei = 1:length(dcp)
    if str2num(dcp{filei}.datapath(end-8:end-1)) < 20220223
        pert3always0deg(filei) = true;
    end
end
if keep_pert3always0deg && sum(pert3always0deg) > 0
    warning('Data set includes at least 1 file where the pertrubation direction during pattern velocity only was always towards 0 deg first, regardless of the direction of pattern velocity')
else
    dcp = dcp(~pert3always0deg);
end
    
%% Get data from initiateCoh experiments
initCohPert = initiateCohPertObj(dcp{1}.sname,dcp{1}.datapath);
initCohPert = initiateCohPertTrials(initCohPert,trialList);

init.t = initCohPert.eye_t';
init.conditions.directions = nan(max(trialList)*length(dcp),1);
init.conditions.speeds = nan(max(trialList)*length(dcp),1);
init.conditions.coh = nan(max(trialList)*length(dcp),1);
init.conditions.pert = nan(max(trialList)*length(dcp),1);
init.eye.hvel = nan(length(init.t),max(trialList)*length(dcp));
init.eye.vvel = nan(length(init.t),max(trialList)*length(dcp));
ind = 1;
for filei = 1:length(dcp)
    
    initCohPert = initiateCohPertObj(dcp{filei}.sname,dcp{filei}.datapath);
    initCohPert = initiateCohPertTrials(initCohPert,trialList);
    
    if ~isempty(initCohPert.trialtype)
        init.conditions.directions(ind:ind+length(initCohPert.directions)-1) = initCohPert.directions;
        init.conditions.speeds(ind:ind+length(initCohPert.directions)-1) = initCohPert.speeds;
        init.conditions.coh(ind:ind+length(initCohPert.directions)-1) = initCohPert.coh;
        init.conditions.pert(ind:ind+length(initCohPert.directions)-1) = initCohPert.perturbations;
        
        init.eye.hvel(:,ind:ind+length(initCohPert.directions)-1) = vertcat(initCohPert.eye(:).hvel)';
        init.eye.vvel(:,ind:ind+length(initCohPert.directions)-1) = vertcat(initCohPert.eye(:).vvel)';
        
        ind = ind+length(initCohPert.directions);
    end
end
init.conditions.directions = init.conditions.directions(1:ind-1);
init.conditions.speeds = init.conditions.speeds(1:ind-1);
init.conditions.coh = init.conditions.coh(1:ind-1);
init.conditions.pert = init.conditions.pert(1:ind-1);
init.eye.hvel = init.eye.hvel(:,1:ind-1);
init.eye.vvel = init.eye.vvel(:,1:ind-1);
init.eye.speed = sqrt(init.eye.hvel.^2 + init.eye.vvel.^2);

%% Find conditional means, marginalize direction and perturbations
speeds = unique(init.conditions.speeds);
cohs = unique(init.conditions.coh);
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        acceptvec = init.conditions.speeds == speeds(si) & init.conditions.coh == cohs(ci) & ...
            ismember(init.conditions.directions,directions);
        init.eye.mean(:,si,ci) = nanmean(init.eye.speed(:,acceptvec),2);
        init.eye.ste(:,si,ci) = nanstd(init.eye.speed(:,acceptvec),[],2)./sqrt(sum(acceptvec));
    end
end


%% Find conditional means, marginalize direction
perturbations = unique(init.conditions.pert);
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        for pi = 1:length(perturbations)
            acceptvec = init.conditions.speeds == speeds(si) & init.conditions.coh == cohs(ci) & init.conditions.pert == perturbations(pi) & ...
                ismember(init.conditions.directions,directions);
            init.eye.pert.mean(:,si,ci,pi) = nanmean(init.eye.speed(:,acceptvec),2);
            init.eye.pert.ste(:,si,ci,pi) = nanstd(init.eye.speed(:,acceptvec),[],2)./sqrt(sum(acceptvec));
            init.eye.pert.N(si,ci,pi) = sum(acceptvec);
            
            % Find perturbation response amplitude
            if perturbations(pi) == 0
                init.eye.pert.t(si,ci,pi) = NaN;
                init.eye.pert.res(si,ci,pi) = 0;
            else

                if perturbations(pi) == 3
                    init.eye.pert.t(si,ci,pi) = 50;
                elseif perturbations(pi) == 7
                    init.eye.pert.t(si,ci,pi) = 600;
                end
                diffTemp = init.eye.pert.mean(init.t >= init.eye.pert.t(si,ci,pi) & ...
                    init.t <= init.eye.pert.t(si,ci,pi)+pertWin,si,ci,pi) - ...
                    init.eye.pert.mean(init.t >= init.eye.pert.t(si,ci,pi) & ...
                    init.t <= init.eye.pert.t(si,ci,pi)+pertWin,si,ci,perturbations == 0);
                steTemp = sqrt(init.eye.pert.ste(init.t >= init.eye.pert.t(si,ci,pi) & ...
                    init.t <= init.eye.pert.t(si,ci,pi)+pertWin,si,ci,pi).^2 + ...
                    init.eye.pert.ste(init.t >= init.eye.pert.t(si,ci,pi) & ...
                    init.t <= init.eye.pert.t(si,ci,pi)+pertWin,si,ci,perturbations == 0).^2);
                [maxRes, maxInd] = max(diffTemp);
                [minRes, minInd] = min(diffTemp);
                init.eye.pert.res(si,ci,pi) = maxRes - minRes;
                init.eye.pert.resSTE(si,ci,pi) = sqrt(steTemp(maxInd).^2 + steTemp(minInd).^2);
                
                if init.eye.pert.t(si,ci,pi) > pertWin
%                     diffTemp2 = init.eye.pert.mean(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
%                         init.t <= init.eye.pert.t(si,ci,pi),si,ci,pi) - ...
%                         init.eye.pert.mean(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
%                         init.t <= init.eye.pert.t(si,ci,pi),si,ci,perturbations == 0);
%                     steTemp2 = sqrt(init.eye.pert.ste(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
%                         init.t <= init.eye.pert.t(si,ci,pi),si,ci,pi).^2 + ...
%                         init.eye.pert.ste(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
%                         init.t <= init.eye.pert.t(si,ci,pi),si,ci,perturbations == 0).^2);
%                     [maxRes, maxInd] = max(diffTemp2);
%                     [minRes, minInd] = min(diffTemp2);
%                     init.eye.pert.resControl(si,ci,pi) = maxRes - minRes;
%                     init.eye.pert.resControlSTE(si,ci,pi) = sqrt(steTemp2(maxInd).^2 + steTemp2(minInd).^2);
                    
                    acceptvec = init.conditions.pert == perturbations(pi) & ...
                        ismember(init.conditions.directions,directions);
                    controlMean = nanmean(init.eye.speed(:,acceptvec),2);
                    controlSTE = nanstd(init.eye.speed(:,acceptvec),[],2)/sqrt(sum(acceptvec));
                    acceptvec = init.conditions.pert == 0 & ...
                        ismember(init.conditions.directions,directions);
                    controlMean2 = nanmean(init.eye.speed(:,acceptvec),2);
                    controlSTE2 = nanstd(init.eye.speed(:,acceptvec),[],2)/sqrt(sum(acceptvec));
                    diffTemp2 = controlMean(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi)) - ...
                        controlMean2(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi));
                    steTemp2 = sqrt(controlSTE(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi)).^2 + ...
                        controlSTE2(init.t >= init.eye.pert.t(si,ci,pi) - pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi)).^2);
                    [maxRes, maxInd] = max(diffTemp2);
                    [minRes, minInd] = min(diffTemp2);
                    init.eye.pert.resControl(si,ci,pi) = maxRes - minRes;
                    init.eye.pert.resControlSTE(si,ci,pi) = sqrt(steTemp2(maxInd).^2 + steTemp2(minInd).^2);
                    
                else
                    acceptvec = init.conditions.pert == perturbations(pi) & ...
                        ismember(init.conditions.directions,directions);
                    controlMean = nanmean(init.eye.speed(:,acceptvec),2);
                    controlSTE = nanstd(init.eye.speed(:,acceptvec),[],2)/sqrt(sum(acceptvec));
                    acceptvec = init.conditions.pert == 0 & ...
                        ismember(init.conditions.directions,directions);
                    controlMean2 = nanmean(init.eye.speed(:,acceptvec),2);
                    controlSTE2 = nanstd(init.eye.speed(:,acceptvec),[],2)/sqrt(sum(acceptvec));
                    diffTemp2 = controlMean(init.t >= init.eye.pert.t(si,ci,pi) + 2*pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi) + 3*pertWin) - ...
                        controlMean2(init.t >= init.eye.pert.t(si,ci,pi) + 2*pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi) + 3*pertWin);
                    steTemp2 = sqrt(controlSTE(init.t >= init.eye.pert.t(si,ci,pi) + 2*pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi) + 3*pertWin).^2 + ...
                        controlSTE2(init.t >= init.eye.pert.t(si,ci,pi) + 2*pertWin & ...
                        init.t <= init.eye.pert.t(si,ci,pi) + 3*pertWin).^2);
                    [maxRes, maxInd] = max(diffTemp2);
                    [minRes, minInd] = min(diffTemp2);
                    init.eye.pert.resControl(si,ci,pi) = maxRes - minRes;
                    init.eye.pert.resControlSTE(si,ci,pi) = sqrt(steTemp2(maxInd).^2 + steTemp2(minInd).^2);                    
                end
            end
            
        end
    end
end

%% Find initiation response
fitvec = [];
testvec = [];
for ci = 1:length(cohs)
    stemp = [];
    for si = 1:length(speeds)
        acceptvec = init.conditions.speeds == speeds(si) & init.conditions.coh == cohs(ci) & ...
            ismember(init.conditions.directions,directions);
        init.eye.init{si,ci} = nanmean(init.eye.speed( init.t>=win(1) & init.t<=win(2),acceptvec),1);
        
        if outlierReject
            for iter = 1:2
                mtemp = nanmean(init.eye.init{si,ci});
                stdtemp = nanstd(init.eye.init{si,ci});
                init.eye.init{si,ci} = init.eye.init{si,ci}(...
                    init.eye.init{si,ci} < mtemp + 1.96*stdtemp & ...
                    init.eye.init{si,ci} > mtemp - 1.96*stdtemp);
            end
        end
        stemp = [stemp; speeds(si)*ones(numel(init.eye.init{si,ci}),1)];
    end
    
    htemp = horzcat(init.eye.init{:,ci})';
    fitvecTemp = false(size(htemp));
    fitInds = randsample(length(htemp),ceil(length(htemp)*0.8));
    fitvecTemp(fitInds) = true;
    [gain(ci).B, gain(ci).BINT, ~, ~, gain(ci).STATS] = ...
        regress(htemp(fitvecTemp),[stemp(fitvecTemp) ones(size(stemp(fitvecTemp)))]);
    
    fitvec = [fitvec; fitvecTemp];
    testvec = [testvec; ~fitvecTemp];
end
gain(1).win = win;

%% Compare regression models of pursuit initiation
res = [];
ss = [];
cs = [];
g = [];
off = [];
fitvec = logical(fitvec);
testvec = logical(testvec);
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        res = [res; init.eye.init{si,ci}'];
        ss = [ss; speeds(si)*ones(size(init.eye.init{si,ci}'))];
        cs = [cs; cohs(ci)*ones(size(init.eye.init{si,ci}'))];
        g = [g; gain(ci).B(1)*ones(size(init.eye.init{si,ci}'))];
        off = [off; gain(ci).B(2)*ones(size(init.eye.init{si,ci}'))];
    end
end
% fitvec = false(size(res));
% fitInds = randsample(length(res),ceil(length(res)*0.8));
% fitvec(fitInds) = true;
% testvec = ~fitvec;
[full.B, full.BINT, ~,~, full.STATS] = regress(res(fitvec),[ss(fitvec) cs(fitvec) ss(fitvec).*cs(fitvec) ones(size(ss(fitvec)))]);
[lin.B, lin.BINT, ~,~, lin.STATS] = regress(res(fitvec),[ss(fitvec) cs(fitvec) ones(size(ss(fitvec)))]);
[interaction.B, interaction.BINT, ~,~, interaction.STATS] = regress(res(fitvec),[ss(fitvec).*cs(fitvec) ones(size(ss(fitvec)))]);

full.rmse = sqrt(mean( (res(testvec) - (full.B'*[ss(testvec) cs(testvec) ss(testvec).*cs(testvec) ones(size(ss(testvec)))]')').^2 ));
lin.rmse = sqrt(mean( (res(testvec) - (lin.B'*[ss(testvec) cs(testvec) ones(size(ss(testvec)))]')').^2 ));
interaction.rmse = sqrt(mean( (res(testvec) - (interaction.B'*[ss(testvec).*cs(testvec) ones(size(ss(testvec)))]')').^2 ));
gain(1).rmse = sqrt(mean( (res(testvec) - (g(testvec).*ss(testvec) + off(testvec))).^2 ));

%% Plotting

%% Mean initiation response
figure;
speedColors = colormap('lines');
cohColors = 1-repmat(cohs/100,[1,3]);
close(gcf)
figure('Name','Mean initCoh response')
tempMax = max([length(cohs) length(speeds)]);
for ci = 1:length(cohs)
    subplot(2,tempMax,ci)
    for si = 1:length(speeds)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = speedColors(si,:);
        myPatch(init.t(:),init.eye.mean(:,si,ci),init.eye.ste(:,si,ci),'patchProperties',patchProps);
        hold on
        plot(init.t,init.eye.mean(:,si,ci),'Color',speedColors(si,:),'LineWidth',2)
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Eye speeds (deg/s)')
    title(['Target coherence = ' num2str(cohs(ci)) ' (deg/s)']);
end

for si = 1:length(speeds)
    subplot(2,tempMax,tempMax+si)
    for ci = 1:length(cohs)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = cohColors(ci,:);
        myPatch(init.t(:),init.eye.mean(:,si,ci),init.eye.ste(:,si,ci),'patchProperties',patchProps);
        hold on
        plot(init.t,init.eye.mean(:,si,ci),'Color',cohColors(ci,:),'LineWidth',2)
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Eye speeds (deg/s)')
    title(['Target speed = ' num2str(speeds(si)) ' (deg/s)']);
end
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

%% Initiation speed
figure('Name','Initiation speed')
xs = linspace(min(speeds),max(speeds),100);
if plotSamps.On
    for ci = 1:length(cohs)
        for si = 1:length(speeds)
            %         plot(speeds(si),randsample(init.eye.init{si,ci},50),'o','Color',colors(ci,:))
            if plotSamps.n < size(init.eye.init{si,ci},2)
                plot(speeds(si)+randn(plotSamps.n,1)/10,init.eye.init{si,ci}(randsample(size(init.eye.init{si,ci},2),plotSamps.n,plotSamps.replacement)),...
                    'o','Color',cohColors(ci,:),'MarkerSize',4,'MarkerFaceColor',cohColors(ci,:))
            else
                plot(speeds(si)+randn(size(init.eye.init{si,ci},2),1)/10,init.eye.init{si,ci},...
                    'o','Color',cohColors(ci,:),'MarkerSize',4,'MarkerFaceColor',cohColors(ci,:))
            end
            hold on
        end
    end
end
for ci = 1:length(cohs)
    for si = 1:length(speeds)
        if plotSamps.On
            plot(speeds(si),nanmean(init.eye.init{si,ci}),...
                'o','Color',cohColors(ci,:),'MarkerSize',8,'MarkerFaceColor',[1 1 1])
        else
            plot(speeds(si),nanmean(init.eye.init{si,ci}),...
                'o','Color',cohColors(ci,:),'MarkerSize',8,'MarkerFaceColor',cohColors(ci,:))
        end
        hold on
    end
    plot(xs,gain(ci).B(1)*xs + gain(ci).B(2),'Color',cohColors(ci,:))
end
plotUnity;
axis square
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])
xlabel('Target speed (deg/s)')
ylabel('Initiation speed (deg/s)')

%% Perturbation response
h = figure;
speedColors = colormap('lines');
set(h,'Position', [347 292 962 937]);
ind = 0;
pertsTemp = perturbations(perturbations ~= 0);
for si = 1:length(speeds)
    for pi = 1:length(pertsTemp)
        ind = ind+1;
        subplot(length(speeds),length(pertsTemp),ind)
        for ci = 1:length(cohs)
            meanCond = init.eye.pert.mean(:,si,ci,perturbations == pertsTemp(pi));
            meanControl = init.eye.pert.mean(:,si,ci,perturbations == 0);
            plot(initCohPert.eye_t,meanCond-meanControl,'Color',cohColors(ci,:));
            hold on
        end
        
        axis tight
        ylim([-speeds(si) speeds(si)]*0.2)
        
        plotHorizontal(0);
        plotVertical(init.eye.pert.t(si,1,perturbations == pertsTemp(pi)));
        
        xlabel('Time from motion onset (ms)')
        ylabel('Perturbation response (deg/s)')
    end
end
ax = axis;
text(ax(2)*0.05,ax(4)*0.95,['Dirs: ' num2str(directions)])

figure('Name','Perturabation response')
pertInd = 0;
for pi = 1:length(perturbations)
    if perturbations(pi) > 0
        pertInd = pertInd+1;
        subplot(1,sum(perturbations>0),pertInd)
        for ci = 1:length(cohs)
            % Gain > 1, therefore one-sided test for significance w/ alpha = 0.05 therefore z = 1.64 (95% chance the true mean is above lower error bar)
            errorbar(speeds',init.eye.pert.res(:,ci,pi)'-init.eye.pert.resControl(:,ci,pi)',...
                sqrt(init.eye.pert.resSTE(:,ci,pi)'.^2 + init.eye.pert.resControlSTE(:,ci,pi)'.^2)*1.64,[],...
                'o-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
            hold on
        end
        xlabel('Target Speed (deg/s)')
        ylabel('Perturbation response')
        title(['Perturbation time = ' num2str(init.eye.pert.t(1,1,pi))])
        axis tight
        ax(pertInd,:) = axis;
    end
end
pertInd = 0;
for pi = 1:length(perturbations)
    if perturbations(pi) > 0
        pertInd = pertInd+1;
        subplot(1,sum(perturbations>0),pertInd)
        axis([0.8*min(speeds) 1.2*max(speeds) min(ax(:,3)) max(ax(:,4))])
        plotHorizontal(0);
    end
end


figure('Name','Feedforward gain estimate')
pertInd = 0;
for pi = 1:length(perturbations)
    if perturbations(pi) > 0
        pertInd = pertInd+1;
        subplot(1,sum(perturbations>0),pertInd)
        for si = 1:length(speeds)
            % Gain > 1, therefore one-sided test for significance w/ alpha = 0.05 therefore z = 1.64 (95% chance the true mean is above lower errorbar)
            errorbar(cohs',(init.eye.pert.res(si,:,pi)'-init.eye.pert.resControl(si,:,pi)')/(0.4*speeds(si)),...
                sqrt(init.eye.pert.resSTE(si,:,pi)'.^2 + init.eye.pert.resControlSTE(si,:,pi)'.^2)/(0.4*speeds(si))*1.64,[],...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            hold on
        end
        xlabel('Coherence')
        ylabel('Gain')
        title(['Perturbation time = ' num2str(init.eye.pert.t(1,1,pi))])
        axis tight
        ax2(pertInd,:) = axis;
    end
end

pertInd = 0;
for pi = 1:length(perturbations)
    if perturbations(pi) > 0
        pertInd = pertInd+1;
        subplot(1,sum(perturbations>0),pertInd)
        axis([0.8*min(cohs) 1.2*max(cohs) min(ax2(:,3)) max(ax2(:,4))])
        plotHorizontal(0);
    end
end

figure('Name','Feedforward gain estimate')
pertInd = 0;
for pi = 1:length(perturbations)
    if perturbations(pi) > 0
        pertInd = pertInd+1;
        subplot(1,sum(perturbations>0),pertInd)
        for ci = 1:length(cohs)
            % Gain > 1, therefore one-sided test for significance w/ alpha = 0.05 therefore z = 1.64 (95% chance the true mean is above lower errorbar)
            errorbar(speeds',(init.eye.pert.res(:,ci,pi)'-init.eye.pert.resControl(:,ci,pi)')./(0.4*speeds'),...
                sqrt(init.eye.pert.resSTE(:,ci,pi)'.^2 + init.eye.pert.resControlSTE(:,ci,pi)'.^2)./(0.4*speeds')*1.64,[],...
                'o-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
            hold on
        end
        xlabel('Target speed (deg/s)')
        ylabel('Gain')
        title(['Perturbation time = ' num2str(init.eye.pert.t(1,1,pi))])
        axis tight
        ax2(pertInd,:) = axis;
    end
end

pertInd = 0;
for pi = 1:length(perturbations)
    if perturbations(pi) > 0
        pertInd = pertInd+1;
        subplot(1,sum(perturbations>0),pertInd)
        axis([0.8*min(speeds) 1.2*max(speeds) min(ax2(:,3)) max(ax2(:,4))])
        plotHorizontal(0);
    end
end
