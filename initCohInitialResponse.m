function initCohInitialResponse(varargin)
%%
%
%
%
%
%%

%% Defaults
speedColors_default = colormap('lines');
cohColors_default = 1-[20 20 20; 60 60 60; 100 100 100]/100;
dynCohRestriction_default.On = true;
dynCohRestriction_default.file = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20220705.mat';
clusterAnalysis_default.On = true;
clusterAnalysis_default.file = '/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohDynCohComp/neuronTyping20220705.mat';

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'binT',0:100:1400)
addParameter(Parser,'dir',[0 180])
addParameter(Parser,'initWin',[100 200])
addParameter(Parser,'win',[-50 250])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'sampleRate',40)
addParameter(Parser,'speedColors',speedColors_default);
addParameter(Parser,'cohColors',cohColors_default);
addParameter(Parser,'clusterAnalysis',clusterAnalysis_default)
addParameter(Parser,'dynCohRestriction',dynCohRestriction_default)

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
binT = Parser.Results.binT;
dir = Parser.Results.dir;
initWin = Parser.Results.initWin;
win = Parser.Results.win;
chanMap = Parser.Results.chanMap;
speedColors = Parser.Results.speedColors;
cohColors = Parser.Results.cohColors;
clusterAnalysis = Parser.Results.clusterAnalysis;
dynCohRestriction = Parser.Results.dynCohRestriction;

%% Load dcp object file
load(dcpObjectFile);


%% Get mean and covariance of each unit
passCutoff = nan(1000,1);
Rinit = nan(length(win(1):win(2)),200,3,3,2,1000);
countsInit = nan(200,3,3,2,1000);
locations = nan(1000,3);
cellID = nan(1000,100,3);
indx = 1;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])

    % Add probe info
    dcp{filei} = addProbeInfo(dcp{filei});

    % InitCoh data
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{filei}.datapath(end-8:end)])

    if ~isempty(initCoh.R)
        if length(initCoh.unitIndex)
            speeds = unique(initCoh.speeds);
            cohs = unique(initCoh.coh);
            dirs = unique(initCoh.directions);

            passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;


            e = sqrt(vertcat(initCoh.eye(:).hvel).^2 + vertcat(initCoh.eye(:).vvel).^2)';
            eInit = nanmean(e(initCoh.eye_t >= initWin(1) & initCoh.eye_t <= initWin(2),:),1);

            % Find counts during initiation window

            % Get data for each neuron
            for uniti = 1:length(initCoh.preferredDirectionRelative)
                ind = find(dir == initCoh.preferredDirectionRelative(uniti));
                antiind = find(dir ~= initCoh.preferredDirectionRelative(uniti));
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(initCoh,dirs(ind),speeds(si),NaN,cohs(ci));
                        rtemp = initCoh.r(:,condLogical,uniti);
                        Rinit(:,1:sum(condLogical),si,ci,1,indx) = rtemp(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2),:);

                        [~,condLogical] = trialSort(initCoh,dirs(antiind),speeds(si),NaN,cohs(ci));
                        rtemp = initCoh.r(:,condLogical,uniti);
                        Rinit(:,1:sum(condLogical),si,ci,2,indx) = rtemp(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2),:);
                    end
                end

                countsTemp = cohConditionedCounts(initCoh,'dirs',dirs(ind),'win',initWin);
                countsInit(1:size(countsTemp,1),:,:,1,indx) = countsTemp(:,:,:,uniti);
                countsTemp = cohConditionedCounts(initCoh,'dirs',dirs(antiind),'win',initWin);
                countsInit(1:size(countsTemp,1),:,:,2,indx) = countsTemp(:,:,:,uniti);

                for j = 1:length(initCoh.unitIndex)
                    cellID(indx,j,1) = filei;
                    cellID(indx,j,2) = initCoh.unitIndex(uniti);
                    cellID(indx,j,3) = initCoh.unitIndex(j);
                end

                if isempty(initCoh.location)
                    x = NaN;
                    y = NaN;
                    z = NaN;
                elseif length(initCoh.location.x)==24
                    siteIndex = chanMap(initCoh.chansIndex(uniti) == chanMap(:,1),2);
                    x = initCoh.location.x(siteIndex);
                    y = initCoh.location.y(siteIndex);
                    depth = -initCoh.location.depth(siteIndex);
                elseif length(initCoh.location.x) > 1
                    siteIndex = floor(initCoh.chansIndex(uniti)/4)+1;
                    tempIndex = find(~isnan(initCoh.location.x));
                    if siteIndex>length(tempIndex)
                        x = NaN;
                        y = NaN;
                        depth = NaN;
                    else
                        x = initCoh.location.x(tempIndex(siteIndex));
                        y = initCoh.location.y(tempIndex(siteIndex));
                        depth = -initCoh.location.depth(tempIndex(siteIndex));
                    end
                else
                    x = initCoh.location.x;
                    y = initCoh.location.y;
                    depth = -initCoh.location.depth;
                end
                locations(indx,:) = [x,y,depth];

                indx = indx+1;
            end
        end
    end
end

Rinit = Rinit(:,:,:,:,:,1:indx-1);
countsInit = countsInit(:,:,:,:,1:indx-1);
locations = locations(1:indx-1,:);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);

%% Remove data that doesn't pass cutoff
countsInit = countsInit(:,:,:,:,passCutoff);
Rinit = Rinit(:,:,:,:,:,passCutoff);
cellID = cellID(passCutoff,:,:);

%% Remove outlier rates
if dynCohRestriction.On
    load(dynCohRestriction.file,'m','m2')
    countsInit = countsInit(:,:,:,:,m<=150 & m2<=150);
    Rinit = Rinit(:,:,:,:,:,m<=150 & m2<=150);
    locations = locations(m<=150 & m2<=150,:);
    cellID = cellID(m<=150 & m2<=150,:,:);
end

%% Find mean spike counts
mCounts = nanmean(countsInit,1);

%% Simple regression model
[Cohs,Speeds] = meshgrid(cohs,speeds);
X = [Speeds(:) Cohs(:) ones(size(Cohs(:)))];
for celli = 1:size(mCounts,4)
    mtemp = reshape(mCounts(1,:,:,1,celli),[length(speeds)*length(cohs) 1]);
    B(:,celli) = regress(mtemp,X);

    mtemp = reshape(mCounts(1,:,:,2,celli),[length(speeds)*length(cohs) 1]);
    antiB(:,celli) = regress(mtemp,X);
end


%% Simple feedforward model comparison

% Generate simple feedforward model predictionmodelN = 1000;
potentialPrefs = 2.^(linspace(-1,8,1000));
tuningParams.N = 10000;
tuningParams.A = 50;
tuningParams.spref = sort(randsample(potentialPrefs,tuningParams.N,true));
tuningParams.sigma = 1.64;
tuningCurve = @(s)( tuningParams.A * exp( -log2(s./tuningParams.spref).^2/(2*tuningParams.sigma.^2) ) );

for si = 1:length(speeds)
    for ci = 1:length(cohs)
        vectorSum(si,ci) = cohs(ci)*sum( tuningCurve(speeds(si)) .* log2(tuningParams.spref) );
    end
end

% Correlate vector sum with firing rates
for uniti = 1:size(mCounts,5)
    tempRates = squeeze(mCounts(1,:,:,1,uniti));
    rtemp = corrcoef(vectorSum(:),tempRates(:));
    vectorSumCorr(uniti) = rtemp(1,2);
end

for uniti = 1:size(Rinit,6)
    for ti = 1:size(Rinit,1)
        tempRates = squeeze(nanmean(Rinit(ti,:,:,:,1,uniti),2));
        rtemp = corrcoef(vectorSum(:),tempRates(:));
        vectorSumCorrT(uniti,ti) = rtemp(1,2);
    end
end


%% Cluster-based analysis
if clusterAnalysis.On
    load(clusterAnalysis.file,'gainRegression','idx','NumClusters')


    for i = 1:NumClusters
        for si = 1:length(speeds)
            for ci = 1:length(cohs)
                clusterRinit(:,si,ci,i) = nanmean(Rinit(:,:,si,ci,1,idx == i),[2,6])*1000;
                mCountsCluster(si,ci,i) = nanmean(countsInit(:,si,ci,1,idx==i),[1,5]);
            end
        end
        
        
        
        for ti = 1:size(clusterRinit,1)
            tempRates = squeeze(clusterRinit(ti,:,:,i));
            rtemp = corrcoef(vectorSum(:),tempRates(:));
            clusterVectorSumCorrT(i,ti) = rtemp(1,2);
        end
        
    end
    
    for ti = 1:size(clusterRinit,1)
        tempRates = squeeze(sum(clusterRinit(ti,:,:,:).*repmat(permute(gainRegression.B(1:NumClusters),[4,2,3,1]),[1,length(speeds),length(cohs),1]),4));
        rtemp = corrcoef(vectorSum(:),tempRates(:));
        weightedVectorSumCorr(i,ti) = rtemp(1,2);
    end

end

%% Compare to behavioral gain during initiation
if testInitGain
    dcpInitCohPert = load('dcpObjectsPertTemp.mat');
    if ~exist('initGain','var')
        [init,~] = initialCohPertBehavioralAnalysis('ar','dcp',dcpInitCohPert.dcp(1:end),...
            'outlierReject',false,'win',[150 200],'keep_pert3always0deg',false,'directions',[0,180],'pertWin',pertWin);
    end
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            gainEst.mean(si,ci) = nanmean(init.eye.init{si,ci}/speeds(si));
            gainEst.ste(si,ci) = nanstd(init.eye.init{si,ci}/speeds(si))/sqrt(numel(init.eye.init{si,ci}));
        end
    end
end

%% Plot results

%% Mean response of example neuron
exIndex = [67, 186];
cellID2 = squeeze(cellID(:,1,1:2));
listIndex = find(ismember(cellID2, exIndex, 'rows'));
figure('Name',['Cell ' num2str(exIndex)],'Position',[312 294 993 960])
subplot(2,2,1)
for ci = 1:length(cohs)

    plot(speeds,squeeze(mCounts(1,:,ci,1,listIndex))/(initWin(2)-initWin(1))*1000,...
        '-','Color',cohColors(ci,:))
    hold on
    for si = 1:length(speeds)
        errorbar(speeds(si),squeeze(mCounts(1,si,ci,1,listIndex))/(initWin(2)-initWin(1))*1000,...
            squeeze(nanstd(countsInit(:,si,ci,1,listIndex),[],1))/sqrt(sum(~isnan(countsInit(:,si,ci,1,listIndex))))/(initWin(2)-initWin(1))*1000,...
            'o-','Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
    end
end
xlim([min(speeds)-5 max(speeds)+5])
yl(1,:) = ylim;

subplot(2,2,3)
for ci = 1:length(cohs)
    plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
        nanmean(Rinit(:,:,2,ci,1,listIndex),2)*1000,'Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
    hold on
    plot(mean(initWin),squeeze(mCounts(1,2,ci,1,listIndex))/(initWin(2)-initWin(1))*1000,'o','Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
end
xlim([min([win(1),initWin(1)]) max([win(2),initWin(2)])])
yl(2,:) = ylim;

subplot(2,2,2)
for ci = 1:length(cohs)

    plot(speeds,squeeze(mCounts(1,:,ci,2,listIndex))/(initWin(2)-initWin(1))*1000,...
        '-','Color',cohColors(ci,:))
    hold on
    for si = 1:length(speeds)
        errorbar(speeds(si),squeeze(mCounts(1,si,ci,2,listIndex))/(initWin(2)-initWin(1))*1000,...
            squeeze(nanstd(countsInit(:,si,ci,2,listIndex),[],1))/sqrt(sum(~isnan(countsInit(:,si,ci,2,listIndex))))/(initWin(2)-initWin(1))*1000,...
            'o-','Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
    end
end
xlim([min(speeds)-5 max(speeds)+5])
yl(3,:) = ylim;

subplot(2,2,4)
for ci = 1:length(cohs)
    plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
        nanmean(Rinit(:,:,2,ci,2,listIndex),2)*1000,'Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
    hold on
    plot(mean(initWin),squeeze(mCounts(1,2,ci,2,listIndex))/(initWin(2)-initWin(1))*1000,'o','Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
end
xlim([min([win(1),initWin(1)]) max([win(2),initWin(2)])])
yl(4,:) = ylim;

for subploti = 1:4
    subplot(2,2,subploti)
    ylim([min(yl(:,1)) max(yl(:,2))])
    if subploti <= 2
        xlabel('Target speed (deg/s)')
        ylabel('Spikes/s')
    else
        plotVertical(initWin);
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
    end

    if subploti == 1
        text(min(speeds)*1.1,min(yl(:,1))*1.1,'Preferred')
    elseif subploti == 2
        text(min(speeds)*1.1,max(yl(:,2))*0.9,'Anti')
    elseif subploti == 3
        if win(1) > 0
            text(win(1)*1.1,min(yl(:,1))*1.1,'Preferred')
        else
            text(win(1)*0.9,min(yl(:,1))*1.1,'Preferred')
        end
    elseif subploti == 4
        if win(1) > 0
            text(win(1)*1.1,max(yl(:,2))*0.9,'Anti')
        else
            text(win(1)*0.9,max(yl(:,2))*0.9,'Anti')
        end
    end
end

%% Example neurons
exUnits = [...
    67,186;...
    69, 81;...
    70,156;...
    75, 80;...
    57,350];
for uniti = 1:size(exUnits,1)

    exIndex = exUnits(uniti,:);
    cellID2 = squeeze(cellID(:,1,1:2));
    listIndex = find(ismember(cellID2, exIndex, 'rows'));
    figure('Name',['Unit ' dcp{exIndex(1)}.datapath(end-11:end) '.' num2str(exIndex(2))],...
        'Position',[1770 1103 751 500])
    subplot(2,2,1)
    for ci = 1:length(cohs)

        plot(speeds,squeeze(mCounts(1,:,ci,1,listIndex))/(initWin(2)-initWin(1))*1000,...
            '-','Color',cohColors(ci,:))
        hold on
        for si = 1:length(speeds)
            errorbar(speeds(si),squeeze(mCounts(1,si,ci,1,listIndex))/(initWin(2)-initWin(1))*1000,...
                squeeze(nanstd(countsInit(:,si,ci,1,listIndex),[],1))/sqrt(sum(~isnan(countsInit(:,si,ci,1,listIndex))))/(initWin(2)-initWin(1))*1000,...
                'o-','Color',cohColors(ci,:),...
                'MarkerFaceColor',cohColors(ci,:))
        end
    end
    xlim([min(speeds)-5 max(speeds)+5])
    yl(1,:) = ylim;

    subplot(2,2,3)
    for ci = 1:length(cohs)
        plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
            nanmean(Rinit(:,:,2,ci,1,listIndex),2)*1000,'Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
        hold on
        plot(mean(initWin),squeeze(mCounts(1,2,ci,1,listIndex))/(initWin(2)-initWin(1))*1000,'o','Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
    end
    xlim([min([win(1),initWin(1)]) max([win(2),initWin(2)])])
    yl(2,:) = ylim;

    subplot(2,2,2)
    for ci = 1:length(cohs)

        plot(speeds,squeeze(mCounts(1,:,ci,2,listIndex))/(initWin(2)-initWin(1))*1000,...
            '-','Color',cohColors(ci,:))
        hold on
        for si = 1:length(speeds)
            errorbar(speeds(si),squeeze(mCounts(1,si,ci,2,listIndex))/(initWin(2)-initWin(1))*1000,...
                squeeze(nanstd(countsInit(:,si,ci,2,listIndex),[],1))/sqrt(sum(~isnan(countsInit(:,si,ci,2,listIndex))))/(initWin(2)-initWin(1))*1000,...
                'o-','Color',cohColors(ci,:),...
                'MarkerFaceColor',cohColors(ci,:))
        end
    end
    xlim([min(speeds)-5 max(speeds)+5])
    yl(3,:) = ylim;

    subplot(2,2,4)
    for ci = 1:length(cohs)
        plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
            nanmean(Rinit(:,:,2,ci,2,listIndex),2)*1000,'Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
        hold on
        plot(mean(initWin),squeeze(mCounts(1,2,ci,2,listIndex))/(initWin(2)-initWin(1))*1000,'o','Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
    end
    xlim([min([win(1),initWin(1)]) max([win(2),initWin(2)])])
    yl(4,:) = ylim;

    for subploti = 1:4
        subplot(2,2,subploti)
        ylim([min(yl(:,1)) max(yl(:,2))])
        if subploti <= 2
            xlabel('Target speed (deg/s)')
            ylabel('Spikes/s')
        else
            plotVertical(initWin);
            xlabel('Time from motion onset (ms)')
            ylabel('Spikes/s')
        end

        if subploti == 1
            text(min(speeds)*1.1,min(yl(:,1))*1.1,'Preferred')
            text(min(speeds)*1.1,min(yl(:,1))*1.3,['Correlation w/ vector sum: ' num2str(vectorSumCorr(listIndex))])
        elseif subploti == 2
            text(min(speeds)*1.1,max(yl(:,2))*0.9,'Anti')
        elseif subploti == 3
            if win(1) > 0
                text(win(1)*1.1,min(yl(:,1))*1.1,'Preferred')
            else
                text(win(1)*0.9,min(yl(:,1))*1.1,'Preferred')
            end
        elseif subploti == 4
            if win(1) > 0
                text(win(1)*1.1,max(yl(:,2))*0.9,'Anti')
            else
                text(win(1)*0.9,max(yl(:,2))*0.9,'Anti')
            end
        end
    end

end

%% Mean spike count across neurons and trials
figure('Name','Average across neurons','Position',[312 294 993 960])
subplot(2,2,1)
for ci = 1:length(cohs)
    plot(speeds,squeeze(nanmean(countsInit(:,:,ci,1,:),[1,5]))/(initWin(2)-initWin(1))*1000,...
        '-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
    hold on

    for si = 1:length(speeds)
        errorbar(speeds(si),squeeze(nanmean(countsInit(:,si,ci,1,:),[1,5]))/(initWin(2)-initWin(1))*1000,...
            squeeze(nanstd(countsInit(:,si,ci,1,:),[],[1,5]))/sqrt(sum(~isnan(countsInit(:,si,ci,1,:)),[1,5]))/(initWin(2)-initWin(1))*1000,...
            'o-','Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
    end
end
xlim([min(speeds)-5 max(speeds)+5])
yl(1,:) = ylim;

subplot(2,2,3)
for ci = 1:length(cohs)
    plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
        nanmean(Rinit(:,:,2,ci,1,:),[2,6])*1000,'Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
    hold on
    plot(mean(initWin),squeeze(nanmean(countsInit(:,2,ci,1,:),[1,5]))/(initWin(2)-initWin(1))*1000,'o','Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
end
xlim([min([win(1),initWin(1)]) max([win(2),initWin(2)])])
yl(2,:) = ylim;

% Antipreferred data
subplot(2,2,2)
for ci = 1:length(cohs)
    plot(speeds,squeeze(nanmean(countsInit(:,:,ci,2,:),[1,5]))/(initWin(2)-initWin(1))*1000,...
        '-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
    hold on

    for si = 1:length(speeds)
        errorbar(speeds(si),squeeze(nanmean(countsInit(:,si,ci,2,:),[1,5]))/(initWin(2)-initWin(1))*1000,...
            squeeze(nanstd(countsInit(:,si,ci,2,:),[],[1,5]))/sqrt(sum(~isnan(countsInit(:,si,ci,2,:)),[1,5]))/(initWin(2)-initWin(1))*1000,...
            'o-','Color',cohColors(ci,:),...
            'MarkerFaceColor',cohColors(ci,:))
    end
end
xlim([min(speeds)-5 max(speeds)+5])
yl(3,:) = ylim;

subplot(2,2,4)
for ci = 1:length(cohs)
    plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
        nanmean(Rinit(:,:,2,ci,2,:),[2,6])*1000,'Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
    hold on
    plot(mean(initWin),squeeze(nanmean(countsInit(:,2,ci,2,:),[1,5]))/(initWin(2)-initWin(1))*1000,'o','Color',cohColors(ci,:),...
        'MarkerFaceColor',cohColors(ci,:))
end
xlim([min([win(1),initWin(1)]) max([win(2),initWin(2)])])
yl(4,:) = ylim;

for subploti = 1:4
    subplot(2,2,subploti)
    ylim([min(yl(:,1)) max(yl(:,2))])
    if subploti <= 2
        xlabel('Target speed (deg/s)')
        ylabel('Spikes/s')
    else
        plotVertical(initWin);
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
    end

    if subploti == 1
        text(min(speeds)*1.1,min(yl(:,1))*1.1,'Preferred')
    elseif subploti == 2
        text(min(speeds)*1.1,max(yl(:,2))*0.9,'Anti')
    elseif subploti == 3
        if win(1) > 0
            text(win(1)*1.1,min(yl(:,1))*1.1,'Preferred')
        else
            text(win(1)*0.9,min(yl(:,1))*1.1,'Preferred')
        end
    elseif subploti == 4
        if win(1) > 0
            text(win(1)*1.1,max(yl(:,2))*0.9,'Anti')
        else
            text(win(1)*0.9,max(yl(:,2))*0.9,'Anti')
        end
    end
end

%% Sort by cluster
if clusterAnalysis.On
    load(clusterAnalysis.file,'gainRegression','idx','NumClusters')

    figure('Name','Initial response sorted by cluster','Position',[182 902 2339 840])
    ind = 0;
    for si = 1:length(speeds)
        for i = 1:NumClusters
            subplot(length(speeds),NumClusters,i+(si-1)*NumClusters)
            for ci = 1:length(cohs)
                plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
                    clusterRinit(:,si,ci,i),'Color',cohColors(ci,:))
                hold on
            end
            axis tight
            ind = ind+1;
            ax(ind,:) = axis;
            title(['Cluster #' num2str(i)])
            xlabel('Time from motion onset (ms)')
            ylabel('Mean firing rate (spikes/s)')
        end
    end
    for i = 1:ind
        subplot(length(speeds),NumClusters,i)
        axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    end


    figure('Name','Weighted average of clusters','Position',[182 902 800 840])
    for si = 1:length(speeds)
        subplot(length(speeds),1,si)
        for ci = 1:length(cohs)
            plot(initCoh.neuron_t(initCoh.neuron_t >= win(1) & initCoh.neuron_t <= win(2)),...
                sum(clusterRinit(:,si,ci,:).*repmat(permute(gainRegression.B(1:NumClusters).*std(gainRegression.x,[],1)',[4,2,3,1]),[size(clusterRinit,1),1,1,1]),4),...
                'Color',cohColors(ci,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Mean firing rate (spikes/s)')
        title(['Speed = ' num2str(speeds(si)) ' deg/s'])
        axis tight
        ylims(si,:) = ylim;
    end
    for si = 1:length(speeds)
        subplot(length(speeds),1,si)
        ylim([min(ylims(:,1)) max(ylims(:,2))])
    end

    figure('Name','Counts as a funciton of speed','Position',[112 373 2410 420])
    for i = 1:NumClusters
        subplot(1,NumClusters+1,i)
        for ci = 1:length(cohs)
            plot(speeds,mCountsCluster(:,ci,i)*1000/(initWin(2)-initWin(1)),'o-',...
                'Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
            hold on
        end
        axis tight
        ylims2(i,:) = ylim;
    end
    subplot(1,NumClusters+1,NumClusters+1)
    for ci = 1:length(cohs)
        plot(speeds,sum(mCountsCluster(:,ci,:).*repmat(permute(gainRegression.B(1:NumClusters).*std(gainRegression.x,[],1)',[3,2,1]),[size(mCountsCluster,1),1,1]),3)*1000/(initWin(2)-initWin(1)),...
            'o-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
        hold on
    end
    xlabel('Target speed (deg/s)')
    ylabel('Spikes/s')
    title(['Weighted average'])
    for i = 1:NumClusters
        subplot(1,NumClusters+1,i)
        ylim([min(ylims2(:,1)) max(ylims2(:,2))])
        xlabel('Target speed (deg/s)')
        ylabel('Spikes/s')
        title(['Cluster #' num2str(i)])
    end
end

%% Compare behavioral gain and neural estimate of gain (from grand mean)
if testInitGain
    figure;
    subplot(1,2,1)
    for ci = 1:length(cohs)
        errorbar(speeds,gainEst.mean(:,ci),gainEst.ste(:,ci),...
            'o-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
        hold on
    end
    axis square
    xlabel('Target speed (deg/s)')
    ylabel('Behavioral gain')
    
    subplot(1,2,2)
    for ci = 1:length(cohs)
        plot(speeds,squeeze(nanmean(countsInit(:,:,ci,1,:),[1,5]))/(initWin(2)-initWin(1))*1000,...
            '-','Color',cohColors(ci,:),'MarkerFaceColor',cohColors(ci,:))
        hold on
        
        for si = 1:length(speeds)
            errorbar(speeds(si),squeeze(nanmean(countsInit(:,si,ci,1,:),[1,5]))/(initWin(2)-initWin(1))*1000,...
                squeeze(nanstd(countsInit(:,si,ci,1,:),[],[1,5]))/sqrt(sum(~isnan(countsInit(:,si,ci,1,:)),[1,5]))/(initWin(2)-initWin(1))*1000,...
                'o-','Color',cohColors(ci,:),...
                'MarkerFaceColor',cohColors(ci,:))
        end
    end
    axis square
    xlabel('Target speed (deg/s)')
    ylabel('Spikes/s')
end