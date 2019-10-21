
subject = 'ar';
dataShort = '191018b';
data = '20191018b';

datapath = ['/home/seth/Projects/DynamicCoherencePhysiology/' subject '/' data];
plxfile = [subject dataShort '.plx'];
plxpath = [datapath(1:end-1) 'plx'];

dcp = dcpObj(subject,datapath);
dcp = extractSpikingData(dcp,plxfile,plxpath,dcp.datapath);

%% Direction preference
dirPref = dirPrefObj(subject,datapath);         % Builds direction preference object
dirPref = assertSpikesExtracted(dirPref,true);  % Asserts that spiking data has been extracted already
dirPref = unitsIndex(dirPref);                  % Finds the indices to the units
dirPref = dirPrefTrials(dirPref,1:2500);        % Finds dirPref trial data

% Sort by direction and plot rasters
dirPrefRaster(dirPref,0:45:315,1);

% Plot grand average PSTH across all trials
figure;
r = dirPref.calcRates(50);                      % Calculates rates based with 50 ms boxcar
plot(-100:1600,mean(r*1000,2));
hold on
plotVertical(0);
xlabel('Time since motion onset (ms)')
ylabel('Sp/s')
mymakeaxis(gca);

%% Speed preference
speedPref = speedPrefObj(subject,datapath);
speedPref = assertSpikesExtracted(speedPref,true);
speedPref = unitsIndex(speedPref);
speedPref = speedPrefTrials(speedPref,1:2500);

% Sort by speed and plot rasters
speedPrefRaster(speedPref,unique(speedPref.directions),unique(speedPref.speeds),1)

% Plot grand average PSTH across all trials
% figure;
% r = speedPref.calcRates(50);                      % Calculates rates based with 50 ms boxcar
% plot(-100:1600,mean(r*1000,2));
% hold on
% plotVertical(0);
% xlabel('Time since motion onset (ms)')
% ylabel('Sp/s')
% mymakeaxis(gca);

%% Dynamic Coherence
dynamicCoh = dynamicCohObj(subject,datapath);
dynamicCoh = assertSpikesExtracted(dynamicCoh,true);
dynamicCoh = unitsIndex(dynamicCoh);
dynamicCoh = dynamicCohTrials(dynamicCoh,1:2500);

% Plot behavioral data
dynamicCohMeanEyeSeq(dynamicCoh,[0 180]);

% Sort by speed and plot rasters
[R,Rste] = dynamicCohSeqConditionedRates(dynamicCoh);