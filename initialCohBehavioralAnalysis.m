function [init, gain] = initialCohBehavioralAnalysis(sname,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'sname')
addParameter(Parser,'dcpObjectsFile',[])
addParameter(Parser,'trialList',1:5000)
addParameter(Parser,'win',[150 200])

parse(Parser,sname,varargin{:})

sname = Parser.Results.sname;
dcpObjectsFile = Parser.Results.dcpObjectsFile;
trialList = Parser.Results.trialList;
win = Parser.Results.win;

%% Build dcp objects
if isempty(dcpObjectsFile)
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
else
    load(dcpObjectsFile,'dcp')
end

%% Get data from initiateCoh experiments
initCoh = initiateCohObj(dcp{1}.sname,dcp{1}.datapath);
initCoh = initiateCohTrials(initCoh,trialList);

init.t = initCoh.eye_t';
init.conditions.directions = nan(max(trialList)*length(dcp),1);
init.conditions.speeds = nan(max(trialList)*length(dcp),1);
init.conditions.coh = nan(max(trialList)*length(dcp),1);
init.eye.hvel = nan(length(init.t),max(trialList)*length(dcp));
init.eye.vvel = nan(length(init.t),max(trialList)*length(dcp));
ind = 1;
for filei = 1:length(dcp)
    
    initCoh = initiateCohObj(dcp{filei}.sname,dcp{filei}.datapath);
    initCoh = initiateCohTrials(initCoh,trialList);
    
    if ~isempty(initCoh.trialtype)
        init.conditions.directions(ind:ind+length(initCoh.directions)-1) = initCoh.directions;
        init.conditions.speeds(ind:ind+length(initCoh.directions)-1) = initCoh.speeds;
        init.conditions.coh(ind:ind+length(initCoh.directions)-1) = initCoh.coh;
        
        init.eye.hvel(:,ind:ind+length(initCoh.directions)-1) = vertcat(initCoh.eye(:).hvel)';
        init.eye.vvel(:,ind:ind+length(initCoh.directions)-1) = vertcat(initCoh.eye(:).vvel)';
        
        ind = ind+length(initCoh.directions);
    end
end
init.conditions.directions = init.conditions.directions(1:ind-1);
init.conditions.speeds = init.conditions.speeds(1:ind-1);
init.conditions.coh = init.conditions.coh(1:ind-1);
init.eye.hvel = init.eye.hvel(:,1:ind-1);
init.eye.vvel = init.eye.vvel(:,1:ind-1);
init.eye.speed = sqrt(init.eye.hvel.^2 + init.eye.vvel.^2);

%% Find conditional means, marginalize direction
speeds = unique(init.conditions.speeds);
cohs = unique(init.conditions.coh);
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        acceptvec = init.conditions.speeds == speeds(si) & init.conditions.coh == cohs(ci);
        init.eye.mean(:,si,ci) = nanmean(init.eye.speed(:,acceptvec),2);
        init.eye.ste(:,si,ci) = nanstd(init.eye.speed(:,acceptvec),[],2)./sqrt(sum(acceptvec));
    end
end

%% Find initiation response
for ci = 1:length(cohs)
    stemp = [];
    for si = 1:length(speeds)
        acceptvec = init.conditions.speeds == speeds(si) & init.conditions.coh == cohs(ci);
        init.eye.init{si,ci} = nanmean(init.eye.speed( init.t>=win(1) & init.t<=win(2),acceptvec),1);
        
        for iter = 1:2
            mtemp = nanmean(init.eye.init{si,ci});
            stdtemp = nanstd(init.eye.init{si,ci});
            init.eye.init{si,ci} = init.eye.init{si,ci}(...
                init.eye.init{si,ci} < mtemp + 1.96*stdtemp & ...
                init.eye.init{si,ci} > mtemp - 1.96*stdtemp);
        end
        stemp = [stemp; speeds(si)*ones(numel(init.eye.init{si,ci}),1)];
    end
    
    [gain(ci).B, gain(ci).BINT, ~, ~, gain(ci).STATS] = ...
        regress(horzcat(init.eye.init{:,ci})',[stemp ones(size(stemp))]);
    
end
gain(1).win = win;

%% Plotting

%% Mean initiation response
colors = colormap('lines');
close(gcf)
figure('Name','Mean initCoh response')
tempMax = max([length(cohs) length(speeds)]);
for ci = 1:length(cohs)
    subplot(2,tempMax,ci)
    for si = 1:length(speeds)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = colors(si,:);
        myPatch(init.t(:),init.eye.mean(:,si,ci),init.eye.ste(:,si,ci),'patchProperties',patchProps);
        hold on
        plot(init.t,init.eye.mean(:,si,ci),'Color',colors(si,:),'LineWidth',2)
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Eye speeds (deg/s)')
    title(['Target coherence = ' num2str(cohs(ci)) ' (deg/s)']);
end

for si = 1:length(speeds)
    subplot(2,tempMax,tempMax+si)
    for ci = 1:length(cohs)
        patchProps.FaceAlpha = 0.3;
        patchProps.FaceColor = colors(ci,:);
        myPatch(init.t(:),init.eye.mean(:,si,ci),init.eye.ste(:,si,ci),'patchProperties',patchProps);
        hold on
        plot(init.t,init.eye.mean(:,si,ci),'Color',colors(ci,:),'LineWidth',2)
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Eye speeds (deg/s)')
    title(['Target speed = ' num2str(speeds(si)) ' (deg/s)']);
end

%% Initiation speed
figure('Name','Initiation speed')
xs = linspace(min(speeds),max(speeds),100);
for ci = 1:length(cohs)
    for si = 1:length(speeds)
%         plot(speeds(si),randsample(init.eye.init{si,ci},50),'o','Color',colors(ci,:))
        hold on
        plot(speeds(si),nanmean(init.eye.init{si,ci}),...
            'o','Color',colors(ci,:),'MarkerSize',6,'MarkerFaceColor',colors(ci,:))
        
        
    end
    plot(xs,gain(ci).B(1)*xs + gain(ci).B(2),'Color',colors(ci,:))
end
plotUnity;
axis square
xlabel('Target speed (deg/s)')
ylabel('Initiation speed (deg/s)')
