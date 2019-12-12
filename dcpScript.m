
subject = 'ar';
dataShort = '191021a';
data = '20191021a';
extractSpikes = true;
extractSpikesKK = false;
addSpike = false;
trials = 1:2500;
unit = 1;
boxCarWidth = 30;
overwrite = true;

datapath = ['/home/seth/Projects/DynamicCoherencePhysiology/' subject '/' data];
if str2num(dataShort(1:end-1)) > 191114
    plxfile = [subject dataShort '.pl2'];
else
    plxfile = [subject dataShort '.plx'];
end
plxpath = [datapath(1:end-1) 'plx'];
kkfile = [subject dataShort '.kwik'];
kkpath = [datapath(1:end-1) 'kk'];

dcp = dcpObj(subject,datapath);
if extractSpikes
    dcp = extractSpikingData(dcp,plxfile,plxpath,dcp.datapath,addSpike);
elseif extractSpikesKK
    dcp = extractKKData(dcp,kkfile,kkpath,plxfile,plxpath,dcp.datapath,addSpike);
end

dcp = tableImport(dcp);

%% Direction preference
dirPref = dirPrefObj(subject,datapath);         % Builds direction preference object
if extractSpikes
    dirPref = assertSpikesExtracted(dirPref,extractSpikes);  % Asserts that spiking data has (not) been extracted already
elseif extractSpikesKK
    dirPref = assertSpikesExtracted(dirPref,extractSpikesKK);  % Asserts that spiking data has (not) been extracted already
end
dirPref = unitsIndex(dirPref);                  % Finds the indices to the units
dirPref = dirPrefTrials(dirPref,trials);        % Finds dirPref trial data

units = [];
for i = 1:length(dirPref.spikeTimes)
    units = [units, dirPref.spikeTimes{i}{2}];
end
units = unique(units);
dirPref.klustaID = units;

% Sort by direction and plot rasters
dirPrefRaster(dirPref,0:45:315,units(unit));

% Plot grand average PSTH across all trials
figure;
r = dirPref.calcRates(boxCarWidth);                      % Calculates rates based with 50 ms boxcar
plot(-100:1600,mean(r(:,:,unit)*1000,2));
hold on
plotVertical(0);
xlabel('Time since motion onset (ms)')
ylabel('Sp/s')
mymakeaxis(gca);

%% Speed preference
speedPref = speedPrefObj(subject,datapath);
if extractSpikes
    speedPref = assertSpikesExtracted(speedPref,extractSpikes);
elseif extractSpikesKK
    speedPref = assertSpikesExtracted(speedPref,extractSpikesKK);
end
speedPref = unitsIndex(speedPref);
speedPref = speedPrefTrials(speedPref,trials);
units = [];
for i = 1:length(speedPref.spikeTimes)
    units = [units, speedPref.spikeTimes{i}{2}];
end
units = unique(units);
speedPref.klustaID = units;

% Sort by speed and plot rasters
speedPrefRaster(speedPref,unique(speedPref.directions),unique(speedPref.speeds),units(unit))

% Plot grand average PSTH across all trials
figure;
r = speedPref.calcRates(50);                      % Calculates rates based with 50 ms boxcar
plot(-100:1600,mean(r(:,:,unit)*1000,2));
hold on
plotVertical(0);
xlabel('Time since motion onset (ms)')
ylabel('Sp/s')
mymakeaxis(gca);

%% Dynamic Coherence
dynamicCoh = dynamicCohObj(subject,datapath);
if extractSpikes
    dynamicCoh = assertSpikesExtracted(dynamicCoh,extractSpikes);
elseif extractSpikesKK
    dynamicCoh = assertSpikesExtracted(dynamicCoh,extractSpikesKK);
end
dynamicCoh = unitsIndex(dynamicCoh);
dynamicCoh = dynamicCohTrials(dynamicCoh,trials);
dirs = [0];

% Plot behavioral data
dynamicCohMeanEyeSeq(dynamicCoh,dirs);
ax = axis;

% Sort by speed and plot rates
t = -100:1600;
[R,Rste] = dynamicCohSeqConditionedRates(dynamicCoh,'width',boxCarWidth,...
    'dirs',dirs,'t',t);

controlInd = 5;
figure;
subplot(2,1,1)
plot(t,R(:,:,unit))
hold on
plotVertical([150 150+0:300:1500]);
xlim(ax(1:2))
subplot(2,1,2)
plot(t,R(:,:,unit) - repmat(R(:,controlInd,unit),[1,size(R,2)]))
hold on
plotVertical([150 150+0:300:1500]);
xlim(ax(1:2))

%% Initiate Coherence
initCoh = initiateCohObj(subject,datapath);
if extractSpikes
    initCoh = assertSpikesExtracted(initCoh,extractSpikes);
elseif extractSpikesKK
    initCoh = assertSpikesExtracted(initCoh,extractSpikesKK);
end
initCoh = unitsIndex(initCoh);
initCoh = initiateCohTrials(initCoh,trials);
dirs = [0];

% Plot behavioral data
initiateCohMeanEye(initCoh,dirs);

% Sort by coherence and plot rates
[R,Rste] = cohConditionedRates(initCoh,'width',boxCarWidth,'dirs',dirs);

controlInd = 5;
figure('Name','Speed/coherence conditioned rates','Position',[342 449 1972 420])
t = -100:1600;
for cohi = 1:length(unique(initCoh.coh))
    subplot(1,length(unique(initCoh.coh)),cohi)
%     for speedi = 1:length(unique(initCoh.speeds))
        plot(t(t<300),R(t<300,:,cohi,unit)*1000)
        %     end
        hold on
        axis tight
        ax(cohi,:) = axis;
end

for cohi = 1:length(unique(initCoh.coh))
    subplot(1,length(unique(initCoh.coh)),cohi)
    axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    plotVertical(0);
    xlabel('Time (ms)')
    ylabel('sp/s')
    mymakeaxis(gca)
end

figure('Name','Speed/coherence conditioned rates 2','Position',[342 449 1972 420])
t = -100:1600;
for speedi = 1:length(unique(initCoh.speeds))
    subplot(1,length(unique(initCoh.speeds)),speedi)
%     for speedi = 1:length(unique(initCoh.speeds))
        plot(t(t<300),squeeze(R(t<300,speedi,:,unit)*1000))
        %     end
        hold on
        axis tight
        ax(speedi,:) = axis;
end

for speedi = 1:length(unique(initCoh.speeds))
    subplot(1,length(unique(initCoh.speeds)),speedi)
    axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
    plotVertical(0);
    xlabel('Time (ms)')
    ylabel('sp/s')
    mymakeaxis(gca)
end