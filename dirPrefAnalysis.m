function dirPrefAnalysis(varargin)
%%
%
%
%
%
%%

%% Defaults
h = figure;
dirColors_default = colormap('lines');
close(h);

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'countWin',[100 200])
addParameter(Parser,'baselineWin',[-50,50])
addParameter(Parser,'win',[-50 250])
addParameter(Parser,'dirs',linspace(0,315,8))
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'sampleRate',40)
addParameter(Parser,'dirColors',dirColors_default)
addParameter(Parser,'boxCarWidth',30)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',10)
addParameter(Parser,'cutWindow',[1 1501]);


parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
countWin = Parser.Results.countWin;
baselineWin = Parser.Results.baselineWin;
win = Parser.Results.win;
dirs = Parser.Results.dirs;
chanMap = Parser.Results.chanMap;
dirColors = Parser.Results.dirColors;
boxCarWidth = Parser.Results.boxCarWidth;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
cutWindow = Parser.Results.cutWindow;

%% Load dcp object file
load(dcpObjectFile);


%% Get mean and covariance of each unit
passCutoff = nan(1000,1);
Rinit = nan(length(win(1):win(2)),200,length(dirs),1000);
counts = nan(200,length(dirs),1000);
baselineCounts = nan(200,length(dirs),1000);
locations = nan(1000,3);
cellID = nan(1000,100,3);

indx = 1;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])

    % Add probe info
    dcp{filei} = addProbeInfo(dcp{filei});
    dcp{filei} = tableImport(dcp{filei});

    % Dir pref data
    dirPref = dirPrefObj(dcp{filei}.sname,dcp{filei}.datapath);
    dirPref = assertSpikesExtracted(dirPref,dcp{filei}.spikesExtracted);
    dirPref = unitsIndex(dirPref);
    dirPref = dirPrefTrials(dirPref,1:5000);
    if ~isempty(dirPref.unitIndex)
        dirPref = dirConditionedRates(dirPref,'width',boxCarWidth);
        dirPref = dirPref.dirRates(boxCarWidth);
        dirPref = dirfindActive(dirPref,rateCutoff,cutWindow);
    end
    dirPref.location = dcp{filei}.location;

    if ~isempty(dirPref.R) 
        if length(dirPref.unitIndex)
            dirs = unique(dirPref.directions);
            
            passCutoff(indx:indx+length(dirPref.unitIndex)-1) = dirPref.passCutoff;

            if length(dirPref.unitIndex) ~= length(dirPref.passCutoff)
                disp(indx)
            end
%             e = sqrt(vertcat(dirPref.eye(:).hvel).^2 + vertcat(dirPref.eye(:).vvel).^2)';
%             eInit = nanmean(e(dirPref.eye_t >= initWin(1) & dirPref.eye_t <= initWin(2),:),1);
            
            % Get data for each neuron
            if checkUnitType && isprop(dcp{filei},'unitTypes')
                unitInd = find(strcmp(dcp{filei}.unitTypes,'good'));
            else
                unitInd = 1:length(dirPref.unitIndex);
            end
            for uniti = unitInd
                
                for di = 1:length(dirs)
                    [~,condLogical] = trialSort(dirPref,dirs(di),NaN,NaN);
                    rtemp = dirPref.r(:,condLogical,uniti);
                    Rinit(:,1:sum(condLogical),di,indx) = rtemp(dirPref.neuron_t >= win(1) & dirPref.neuron_t <= win(2),:);
                end
                
                countsTemp = dirConditionedCounts(dirPref,'dirs',dirs,'win',countWin);
                counts(1:size(countsTemp,1),:,indx) = countsTemp(:,:,uniti);
                
                countsTemp = dirConditionedCounts(dirPref,'dirs',dirs,'win',baselineWin);
                baselineCounts(1:size(countsTemp,1),:,indx) = countsTemp(:,:,uniti);
                
                for j = 1:length(dirPref.unitIndex)
                    cellID(indx,j,1) = filei;
                    cellID(indx,j,2) = dirPref.unitIndex(uniti);
                    cellID(indx,j,3) = dirPref.unitIndex(j);
                end

                if isempty(dirPref.location)
                    x = NaN;
                    y = NaN;
                    depth = NaN;
                elseif length(dirPref.location.x)==24
                    siteIndex = chanMap(dirPref.chansIndex(uniti) == chanMap(:,1),2);
                    x = dirPref.location.x(siteIndex);
                    y = dirPref.location.y(siteIndex);
                    depth = -dirPref.location.depth(siteIndex);
                elseif length(dirPref.location.x) > 1
                    siteIndex = floor(dirPref.chansIndex(uniti)/4)+1;
                    tempIndex = find(~isnan(dirPref.location.x));
                    if siteIndex>length(tempIndex)
                        x = NaN;
                        y = NaN;
                        depth = NaN;
                    else
                        x = dirPref.location.x(tempIndex(siteIndex));
                        y = dirPref.location.y(tempIndex(siteIndex));
                        depth = -dirPref.location.depth(tempIndex(siteIndex));
                    end
                else
                    x = dirPref.location.x;
                    y = dirPref.location.y;
                    depth = -dirPref.location.depth;
                end
                locations(indx,:) = [x,y,depth];

                indx = indx+1;
            end
        end
    end
end

Rinit = Rinit(:,:,:,1:indx-1);
counts = counts(:,:,1:indx-1);
baselineCounts = baselineCounts(:,:,1:indx-1);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);
locations = locations(1:indx-1,:);

%% Remove data that doesn't pass cutoff
counts = counts(:,:,passCutoff);
baselineCounts = baselineCounts(:,:,passCutoff);
Rinit = Rinit(:,:,:,passCutoff);
cellID = cellID(passCutoff,:,:);
locations = locations(passCutoff,:);

locations(locations(:,3) < -2e4) = NaN;

%% Find mean spike counts
mCounts = nanmean(counts,1);

%% direction index
DI = (max(mCounts,[],2) - min(mCounts,[],2))./(max(mCounts,[],2) + min(mCounts,[],2));
DI2 = (max(mCounts,[],2) - min(mCounts,[],2))./(max(mCounts,[],2));

%% Find preferred direction
[~,maxInd] = max(mCounts,[],2);
preferredDirection = dirs(squeeze(maxInd));

% K nearest neighbor distances
locations2 = locations;
locations2(locations(:,1)>1,:) = [locations(locations(:,1)>1,2), locations(locations(:,1)>1,1), locations(locations(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
locations2(:,3) = locations2(:,3)/1000;
Kneighbors = 2;
nnIdx = knnsearch(locations2,locations2,'K',Kneighbors+1);
dirDif = nan(size(nnIdx,1),Kneighbors);
dirDifRand = nan(size(nnIdx,1),Kneighbors);
aveDist = nan(size(nnIdx,1),1);
aveDistRand = nan(size(nnIdx,1),1);
for ni = 1:size(nnIdx,1)
    aveDist(ni) = mean( sqrt( sum((repmat(locations2(ni,:),[Kneighbors,1]) - locations2(nnIdx(ni,2:end),:)).^2,2) ) );
    dirDif(ni,:) = wrapTo180(preferredDirection(ni) - preferredDirection(nnIdx(ni,2:end)));
    randInds = randsample(size(nnIdx,1),Kneighbors);
    dirDifRand(ni,:) = wrapTo180(preferredDirection(ni) - preferredDirection(randInds));
    aveDistRand(ni) = mean( sqrt( sum((repmat(locations2(ni,:),[Kneighbors,1]) - locations2(randInds,:)).^2) ) );
end

distances = pdist2(locations2,locations2);
mask = triu(ones(size(distances)),1);
mask(mask ~= 1) = NaN;
distances = mask.*distances;

%% Direction responses of example neuron
switch dirPref.sname
    case 'ar'
        exIndex = [67, 186];
    case 'fr'
        exIndex = [109, 20];
end
cellID2 = squeeze(cellID(:,1,1:2));
listIndex = find(ismember(cellID2, exIndex, 'rows'));
figure('Name',['Cell ' num2str(exIndex)],'Position',[312 294 993 960])
subplot(1,2,1)
for di = 1:length(dirs)
    plot(win(1):win(2),nanmean(Rinit(:,:,di,listIndex),2)*1000,'Color',dirColors(di,:))
    hold on
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')
title(['Cell ' num2str(exIndex)])
ax = axis;
text((ax(2)-ax(1))*.05+ax(1),(ax(4)-ax(3))*.9+ax(3),...
    [dcp{exIndex(1)}.datapath(end-11:end) '.' num2str(exIndex(2))])

subplot(1,2,2)
thetas = linspace(-pi,pi,60);
polarplot(thetas,...
    nanmean(baselineCounts(:,:,listIndex),[1,2])'*1000/(baselineWin(2)-baselineWin(1))+ones(size(thetas))*1.96*nanstd(baselineCounts(:,:,listIndex),[],[1,2])/sqrt(sum(~isnan(baselineCounts(:,:,listIndex)),[1,2]))*1000/(baselineWin(2)-baselineWin(1)),...
    'k:')
hold on
polarplot(thetas,...
    nanmean(baselineCounts(:,:,listIndex),[1,2])'*1000/(baselineWin(2)-baselineWin(1))-ones(size(thetas))*1.96*nanstd(baselineCounts(:,:,listIndex),[],[1,2])/sqrt(sum(~isnan(baselineCounts(:,:,listIndex)),[1,2]))*1000/(baselineWin(2)-baselineWin(1)),...
    'k:')
for di = 1:length(dirs)
    polarplot(dirs(di)*pi/180,nanmean(counts(:,di,listIndex),1)'*1000/(countWin(2)-countWin(1)),'o','Color',dirColors(di,:),'MarkerFaceColor',dirColors(di,:))
end
countstemp = nanmean(counts(:,:,listIndex),1)*1000/(countWin(2)-countWin(1));
polarplot([dirs; dirs(1)]*pi/180,[countstemp countstemp(1)]','k-')

%% Direction preference and index distribution
figure('Name','Direction index distribution')
subplot(1,2,1)
polarhistogram(preferredDirection*pi/180,linspace(-pi/8,2*pi-pi/8,9))
% xlabel('Preferred direction')
% ylabel('Number of neurons')

subplot(1,2,2)
histogram(DI,linspace(0,1,100));
xlabel('Direction index')
ylabel('Number of neurons')

%% Topography of directions
for di = 1:length(dirs)
    hue = dirs(di)/ 180 / pi;
    saturation = ones(size(hue));
    value = ones(size(hue));
    hsv = cat(3, hue, saturation, value);
    colorWheel(di,:) = hsv2rgb(hsv);
end
colorWheel(colorWheel > 1) = 1;
colorWheel(colorWheel < 0) = 0;

figure;
randScale = 0.08;
vecInds = [1,2,3];
locations2 = locations;
locations2(locations(:,1)>1,:) = [locations(locations(:,1)>1,2), locations(locations(:,1)>1,1), locations(locations(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
locationsRand = locations2 + ...
    [randScale*nanstd(locations2(:,1))*randn(size(locations,1),1), ...
    randScale*nanstd(locations2(:,2))*randn(size(locations,1),1), ...
    0*nanstd(locations2(:,3))*randn(size(locations,1),1)];                  % Add randomness to a-p and m-l locations to make
for uniti = 1:size(locations,1)
    plot3(locationsRand(uniti,1),locationsRand(uniti,2),locationsRand(uniti,3)/1000,...
        'o','Color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',colorWheel(maxInd(1,1,uniti),:,:));
    hold on
end
grid on
axis equal
xlabel('Anterior/postieror (mm)')
ylabel('Medial/lateral (mm)')
zlabel('Depth (mm)')

%% Distribution of distances and direction preferrence differences
figure('Name','Functional vs physical topography','Position',[1153 924 1373 405])
subplot(1,2,1)
histogram(aveDist,linspace(0,max([aveDist; aveDistRand]),50))
hold on
histogram(aveDistRand,linspace(0,max([aveDist; aveDistRand]),50))
xlabel('Average distance (mm)')
ylabel('N')
legend({[num2str(Kneighbors) ' nearest in chamber'],'Random'})

subplot(1,2,2)
histogram(dirDif(:),-202.5:45:202.5)
hold on
histogram(dirDifRand(:),-202.5:45:202.5)
xlabel('Difference in preferred direction (deg)')
ylabel('N')
legend({[num2str(Kneighbors) ' nearest in chamber'],'Random'})

%% Distribution of direction preference differences, neurons w/in X mm
X = logspace(log10(0.05),log10(50),6);
dP = preferredDirection - preferredDirection';
figure('Name','Direction preference difference at different cortical distances',...
    'Position',[1956 97 570 1225])
for i = 1:length(X)
    local = dP(distances <= X(i));
    subplot(length(X),1,i)
    histogram(local,-202.5:45:202.5)
    xlabel('Difference in direction preferrence (deg)')
    ylabel(['Neurons within ' num2str(X(i)*1000) ' um'])
end

%% Plot direction preferences of neurons
for neuroni = 1:size(Rinit,4)
    cellID2 = squeeze(cellID(:,1,1:2));
    figure('Name',['Cell ' num2str(cellID2(neuroni,:))],'Position',[312 294 993 960])
    subplot(1,2,1)
    for di = 1:length(dirs)
        plot(win(1):win(2),nanmean(Rinit(:,:,di,neuroni),2)*1000,'Color',dirColors(di,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    title(['Cell ' num2str(cellID2(neuroni,:))])
    ax = axis;
    text((ax(2)-ax(1))*.05+ax(1),(ax(4)-ax(3))*.9+ax(3),...
        [dcp{cellID2(neuroni,1)}.datapath(end-11:end) '.' num2str(cellID2(neuroni,2))])
    
    subplot(1,2,2)
    thetas = linspace(-pi,pi,60);
    polarplot(thetas,...
        ones(size(thetas))*nanmean(baselineCounts(:,:,neuroni),[1,2])'*1000/(baselineWin(2)-baselineWin(1))+1.96*nanstd(baselineCounts(:,:,neuroni),[],[1,2])/sqrt(sum(~isnan(baselineCounts(:,:,neuroni)),[1,2]))*1000/(baselineWin(2)-baselineWin(1)),...
        'k:')
    hold on
    polarplot(thetas,...
        ones(size(thetas))*nanmean(baselineCounts(:,:,neuroni),[1,2])'*1000/(baselineWin(2)-baselineWin(1))-1.96*nanstd(baselineCounts(:,:,neuroni),[],[1,2])/sqrt(sum(~isnan(baselineCounts(:,:,neuroni)),[1,2]))*1000/(baselineWin(2)-baselineWin(1)),...
        'k:')
    for di = 1:length(dirs)
        polarplot(dirs(di)*pi/180,nanmean(counts(:,di,neuroni),1)'*1000/(countWin(2)-countWin(1)),...
            'o','Color',dirColors(di,:),'MarkerFaceColor',dirColors(di,:))
    end
    countstemp = nanmean(counts(:,:,neuroni))*1000/(countWin(2)-countWin(1));
    polarplot([dirs; dirs(1)]*pi/180,[countstemp countstemp(1)],'k-')
    
    input('Press enter to continue...')
end