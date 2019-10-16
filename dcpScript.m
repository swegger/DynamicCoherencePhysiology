
subject = 'ar';
dataShort = '191015a';
data = '20191015a';

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