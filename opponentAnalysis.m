function opponentAnalysis(varargin)
%%
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'binT',0:100:1400)
addParameter(Parser,'dir',[0 180])
addParameter(Parser,'initWin',[150 200])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'win',-200:200)
addParameter(Parser,'calcCC',false)
addParameter(Parser,'calcRSC',true)
addParameter(Parser,'shuffleN',100)
addParameter(Parser,'sampleRate',40)
addParameter(Parser,'ClusterMethod','densityClust')
addParameter(Parser,'pertWin',250)
addParameter(Parser,'resultsFile','none')
addParameter(Parser,'testInitGain',true)
addParameter(Parser,'dcpInitCohPertFile','dcpObjectsPertTemp.mat')
addParameter(Parser,'dcpDynCohFile',[])
addParameter(Parser,'initSpeed',10)
addParameter(Parser,'includeEarlyInitCohPertTime',false)
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'sequences',[1; 2; 3; 4; 5])
addParameter(Parser,'perturbations',[0; 4; 6; 8])
addParameter(Parser,'NumClusters',8)

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
binT = Parser.Results.binT;
dir = Parser.Results.dir;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
win = Parser.Results.win;
calcCC = Parser.Results.calcCC;
calcRSC = Parser.Results.calcRSC;
shuffleN = Parser.Results.shuffleN;
sampleRate = Parser.Results.sampleRate;
ClusterMethod = Parser.Results.ClusterMethod;
pertWin = Parser.Results.pertWin;
resultsFile = Parser.Results.resultsFile;
testInitGain = Parser.Results.testInitGain;
initSpeed = Parser.Results.initSpeed;
includeEarlyInitCohPertTime = Parser.Results.includeEarlyInitCohPertTime;
dcpInitCohPertFile = Parser.Results.dcpInitCohPertFile;
dcpDynCohFile = Parser.Results.dcpDynCohFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
sequences = Parser.Results.sequences;
perturbations = Parser.Results.perturbations;
NumClusters = Parser.Results.NumClusters;


%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;


%% Load dcp object file
if strcmp(resultsFile,'none')
    load(dcpObjectFile);
    fileExist = false;
else
    fileExist = exist(resultsFile,'file');
    if fileExist
        load(resultsFile)
    else
        error(['File ' resultsFile ' does not exist in path.'])
    end
end

%% Find the gain from behavioral data

if ~exist('gain','var')
    disp('Determining gain from dynamic coherence trials...')
    if isempty(dcpDynCohFile)
        [~,gain] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcp,'directions',[0,180],'pertWin',pertWin);
    else
        dcpDynCoh = load(dcpDynCohFile);
        [~,gain] = dynamicCohBehavioralAnalysis(dcp{1}.sname,'dcp',dcpDynCoh.dcp,'directions',[0,180],'pertWin',pertWin);
    end
end

if testInitGain
    disp('Determining gain from initiate coherence trials...')
    dcpInitCohPert = load(dcpInitCohPertFile);
    if ~exist('initGain','var')
        [init,~] = initialCohPertBehavioralAnalysis(dcp{1}.sname,'dcp',dcpInitCohPert.dcp(1:end),...
            'outlierReject',false,'win',[150 200],'keep_pert3always0deg',false,'directions',[0,180],'pertWin',pertWin);
        initGain = (init.eye.pert.res - init.eye.pert.resControl)./(0.4*repmat(speeds,[1,3,3]));
    end
end

%% Get mean and covariance of each unit

disp('Neural analysis loop...')

passCutoff = nan(1000,1);
Rinit = nan(1701,3,3,1000,2);
Rdyn = nan(1701,5,1000,2);
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
    
    % DynCoh data
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    if ~isempty(initCoh.R) && (~isempty(dynCoh.R) || ~any(isnan(dynCoh.R(:))) && size(dynCoh.R,2) > 1)
        if length(initCoh.unitIndex) == length(dynCoh.unitIndex)
            
            passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff | dynCoh.passCutoff;
            
            
            e = sqrt(vertcat(initCoh.eye(:).hvel).^2 + vertcat(initCoh.eye(:).vvel).^2)';
            eInit = nanmean(e(initCoh.eye_t >= initWin(1) & initCoh.eye_t <= initWin(2),:),1);
                        
            % Get data for each neuron
            for uniti = 1:length(initCoh.preferredDirectionRelative)
                ind = find(dir == initCoh.preferredDirectionRelative(uniti));
                Rinit(:,:,:,indx,1) = initCoh.R(:,:,:,uniti,ind);
                ind = find(dir ~= initCoh.preferredDirectionRelative(uniti));
                Rinit(:,:,:,indx,2) = initCoh.R(:,:,:,uniti,ind);
                
                ind = find(dir == dynCoh.preferredDirectionRelative(uniti));
                Rdyn(:,:,indx,1) = dynCoh.R(:,:,uniti,ind);
                ind = find(dir ~= dynCoh.preferredDirectionRelative(uniti));
                Rdyn(:,:,indx,2) = dynCoh.R(:,:,uniti,ind);
                
                
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

Rinit = Rinit(:,:,:,1:indx-1,:);
Rdyn = Rdyn(:,:,1:indx-1,:);
locations = locations(1:indx-1,:);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff,:);
Rdyn = Rdyn(:,:,passCutoff,:);
locations = locations(passCutoff,:);
cellID = cellID(passCutoff,:,:);

%% Remove outlier rates
m = squeeze(max(Rinit,[],[1,2,3,5]))*1000;
m2 = squeeze(max(Rdyn,[],[1,2,4]))*1000;
Rinit = Rinit(:,:,:,m<=150 & m2<=150,:);
Rdyn = Rdyn(:,:,m<=150 & m2<=150,:);
locations = locations(m<=150 & m2<=150,:);
cellID = cellID(m<=150 & m2<150,:,:);

%% Remove tail
Rdyn = Rdyn(dynCoh.neuron_t<=1350,:,:,:);

%% Mean center
Rinit2 = Rinit;
mRinit2 = nanmean(Rinit2,[1,2,3,5]);
Rinit2 = Rinit2 - mRinit2;

Rdyn2 = Rdyn;
mRdyn2 = nanmean(Rdyn2,[1,2,4]);
Rdyn2 = Rdyn2 - mRdyn2;

%% Reshape
sz = size(permute(Rinit2,[1,2,3,5,4]));
Rinit3 = reshape(permute(Rinit2,[1,2,3,5,4]),[prod(sz(1:4)),prod(sz(5))]);

sz = size(permute(Rdyn2,[1,2,4,3]));
Rdyn3 = reshape(permute(Rdyn2,[1,2,4,3]),[prod(sz(1:3)),prod(sz(end))]);

%% PCA
CR = cov(Rinit3);
[vecs,vals] = eig(CR);
vecs = fliplr(vecs);
vals = flipud(diag(vals));
reconstructionN = 5;

sz = size(Rinit);
RlowD = nan([sz(1:3),reconstructionN,2]);
for si = 1:size(Rinit,2)
    for ci = 1:size(Rinit,3)
        RlowD(:,si,ci,:,1) = permute(vecs(:,1:reconstructionN)'*squeeze(Rinit2(:,si,ci,:,1))',[2,3,4,1]);
        RlowD(:,si,ci,:,2) = permute(vecs(:,1:reconstructionN)'*squeeze(Rinit2(:,si,ci,:,2))',[2,3,4,1]);
    end
end

sz = size(Rdyn);
RDynlowD = nan([sz(1:2),reconstructionN,2]);
for seqi = 1:size(Rdyn,2)
    
    RDynlowD(:,seqi,:,1) = permute(vecs(:,1:reconstructionN)'*squeeze(Rdyn2(:,seqi,:,1))',[2,3,4,1]);
    RDynlowD(:,seqi,:,2) = permute(vecs(:,1:reconstructionN)'*squeeze(Rdyn2(:,seqi,:,2))',[2,3,4,1]);
    
end

%% Simple avaraging
a = Rinit(dynCoh.neuron_t<=1350,:,:,:,:);
b = nan(size(Rdyn));
c = nan(size(Rdyn));
for seqi = 1:size(Rdyn,2)
    t20 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 20);
    if ~isempty(t20)
        b(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:,:) = ...
            Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:,:) - ...
            Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),5,:,:);
    end

    t100 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 100);
    if ~isempty(t100)
        c(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:,:) = ...
            Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:,:) - ...
            Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),5,:,:);
    end
end

figure
subplot(2,2,1)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(c(:,seqi,:,1)),2)*1000,'Color',colors(seqi,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(b(:,seqi,:,1)),2)*1000,'Color',colors(seqi,:))
end
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,1,:,1)-a(:,2,2,:,1)),2)*1000,'Color',initColors(1,:))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,3,:,1)-a(:,2,2,:,1)),2)*1000,'Color',initColors(3,:))
plotHorizontal(0);
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Excess spike/s')

subplot(2,2,2)
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,2,:,1)),2)*1000,'Color',initColors(2,:))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,5,:,1)),2)*1000,'Color',colors(5,:))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

subplot(2,2,3)
for seqi = 1:5
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(c(:,seqi,:,2)),2)*1000,'Color',colors(seqi,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(b(:,seqi,:,2)),2)*1000,'Color',colors(seqi,:))
end
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,1,:,2)-a(:,2,2,:,2)),2)*1000,'Color',initColors(1,:))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,3,:,2)-a(:,2,2,:,2)),2)*1000,'Color',initColors(3,:))
plotHorizontal(0);
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Excess spike/s')

subplot(2,2,4)
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,2,2,:,2)),2)*1000,'Color',initColors(2,:))
hold on
plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,5,:,2)),2)*1000,'Color',colors(5,:))
plotVertical([450 750 1050]);
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

figure('Name','Grand average','Position',[481 159 862 982])
minMax = [Inf,0];
for speedi = 1:length(speeds)
    subplot(length(speeds),2,1+(speedi-1)*2)
    for cohi = 1:length(cohs)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speedi,cohi,:,1)),2)*1000,...
            'Color',[speedColors(speedi,:), cohs(cohi)/100])
        hold on
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speedi,cohi,:,2)),2)*1000,...
            '--','Color',[speedColors(speedi,:), cohs(cohi)/100])
    end
    temp = ylim;
    if minMax(1) > temp(1)
        minMax(1) = temp(1);
    end
    if minMax(2) < temp(2)
        minMax(2) = temp(2);
    end
end
for speedi = 1:length(speeds)
    subplot(length(speeds),2,1+(speedi-1)*2)
    ylim(minMax)
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
end

speedi = find(speeds == 10);
subplot(length(speeds),2,2 + (speedi-1)*2)
for seqi = 1:length(sequences)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,seqi,:,1)),2)*1000,'Color',colors(seqi,:))
    hold on
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,seqi,:,2)),2)*1000,'--','Color',colors(seqi,:))
end
ylim(minMax)
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')


figure('Name','Grand average opponent response','Position',[481 159 862 982])
minMax = [Inf,0];
for speedi = 1:length(speeds)
    subplot(length(speeds),2,1+(speedi-1)*2)
    for cohi = 1:length(cohs)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(-diff(a(:,speedi,cohi,:,:),1,5)),2)*1000,...
            'Color',[speedColors(speedi,:), cohs(cohi)/100])
        hold on
    end
    temp = ylim;
    if minMax(1) > temp(1)
        minMax(1) = temp(1);
    end
    if minMax(2) < temp(2)
        minMax(2) = temp(2);
    end
end
for speedi = 1:length(speeds)
    subplot(length(speeds),2,1+(speedi-1)*2)
    ylim(minMax)
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
end

speedi = find(speeds == 10);
subplot(length(speeds),2,2 + (speedi-1)*2)
for seqi = 1:length(sequences)
    plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(-diff(Rdyn(:,seqi,:,:),1,4)),2)*1000,'Color',colors(seqi,:))
    hold on
end
ylim(minMax)
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

gvmrh = figure('Name',['Grand average opponent response']);
dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,:,1) - Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,:,2),[1,3]))*1000;
for seqi = 1:length(sequences)
    plot(gain(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
    hold on
end
dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,:,1) - Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,:,2),[1,3]))*1000;
for seqi = 1:length(sequences)
    plot(gain(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
end
dynRatesTemp = squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,:,1) - Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,:,2),[1,3]))*1000;
for seqi = 1:length(sequences)
    plot(gain(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:))
end
if testInitGain
    for si = 1:length(speeds)
        initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:,1) - a(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,si,:,:,2),[1,4]))*1000;
        plot(squeeze(initGain(si,:,3)),initRatesTemp,...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
    end
    if includeEarlyInitCohPertTime
        for si = 1:length(speeds)
            initRatesTemp = squeeze(nanmean(a(dynCoh.neuron_t >= 50 & dynCoh.neuron_t <= 150,si,:,:,1) - a(dynCoh.neuron_t >= 50 & dynCoh.neuron_t <= 150,si,:,:,2),[1,4]))*1000;
            plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
        end
    end
    
end
xlabel('Behavioral gain (unitless)')
ylabel('Spikes/s')

%% PCA plotting

figure
for ri = 1:reconstructionN
    for si = 1:size(RlowD,2)
        subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
        for ci = 1:size(RlowD,3)
            plot(initCoh.neuron_t,RlowD(:,si,ci,ri,1),'Color',initColors(ci,:))
            hold on
            plot(initCoh.neuron_t,RlowD(:,si,ci,ri,2),'--','Color',initColors(ci,:))
        end
        templims(si,:) = axis;
    end

    for si = 1:size(RlowD,2)
        subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
        ylim([min(templims(:,3)) max(templims(:,4))]);
    end
        
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
end

figure
for si = 1:size(RlowD,2)
    for ci = 1:size(RlowD,3)
        plot3(RlowD(:,si,ci,1,1),RlowD(:,si,ci,2,1),RlowD(:,si,ci,3,1),'-','Color',[speedColors(si,:) cohs(ci)/100])
        hold on
        plot3(RlowD(:,si,ci,1,2),RlowD(:,si,ci,2,2),RlowD(:,si,ci,3,2),'--','Color',[speedColors(si,:) cohs(ci)/100])
        scatter3(RlowD(1,si,ci,1,1),RlowD(1,si,ci,2,1),RlowD(1,si,ci,3,1),100,speedColors(si,:),...
            'filled','o','MarkerEdgeAlpha',cohs(ci)/100,'MarkerFaceAlpha',cohs(ci)/100)
        scatter3(RlowD(end,si,ci,1,1),RlowD(end,si,ci,2,1),RlowD(end,si,ci,3,1),100,speedColors(si,:),...
            'filled','s','MarkerEdgeAlpha',cohs(ci)/100,'MarkerFaceAlpha',cohs(ci)/100)
        scatter3(RlowD(1,si,ci,1,2),RlowD(1,si,ci,2,2),RlowD(1,si,ci,3,2),100,speedColors(si,:),...
            'filled','o','MarkerEdgeAlpha',cohs(ci)/100,'MarkerFaceAlpha',cohs(ci)/100)
        scatter3(RlowD(end,si,ci,1,2),RlowD(end,si,ci,2,2),RlowD(end,si,ci,3,2),100,speedColors(si,:),...
            'filled','s','MarkerEdgeAlpha',cohs(ci)/100,'MarkerFaceAlpha',cohs(ci)/100)
    end
end
grid on
xlabel('PC_1')
ylabel('PC_2')
zlabel('PC_3')
view([-75 60])