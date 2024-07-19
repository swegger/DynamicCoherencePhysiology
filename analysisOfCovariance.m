function analysisOfCovariance(subjects,varargin)
%%
%
%
%
%
%%

%% Defaults
plotOpts_default.On = false;

calcCov_default.On = true;
calcCov_default.binT = 0:100:1000;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subjects')
addParameter(Parser,'dcpObjectFile',[])
addParameter(Parser,'sourceDirectory',[])
addParameter(Parser,'initCohCollate',true)
addParameter(Parser,'dynCohCollate',false)
addParameter(Parser,'FunctionalTopologyFile','/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat')
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'chanMap',{[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24],...
    [23 1; 22 2; 21 3; 20 4; 19 5; 18 6; 17 7; 16 8; 15 9; 14 10; 13 11; 12 12; 11 13; 10 14; 9 15; 8 16; 7 17; 6 18; 5 19; 4 20; 3 21; 2 22; 1 23; 0 24]})
addParameter(Parser,'calcCov',calcCov_default)
addParameter(Parser,'pertWin',250)
addParameter(Parser,'dcpInitCohPertFile','dcpObjectsPertTemp.mat')
addParameter(Parser,'initSpeed',10)
addParameter(Parser,'includeEarlyInitCohPertTime',true)
addParameter(Parser,'dcpDynCohFile',{[]})
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'sequences',[1; 2; 3; 4; 5])
addParameter(Parser,'NumClusters',8)
addParameter(Parser,'checkUnitType',false)
addParameter(Parser,'rateCutoff',NaN)
addParameter(Parser,'winT',[0 1000])
addParameter(Parser,'plotOpts',plotOpts_default)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,subjects,varargin{:})

subjects = Parser.Results.subjects;
dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
initCohCollate = Parser.Results.initCohCollate;
dynCohCollate = Parser.Results.dynCohCollate;
FunctionalTopologyFile = Parser.Results.FunctionalTopologyFile;
directions = Parser.Results.directions;
chanMap = Parser.Results.chanMap;
calcCov = Parser.Results.calcCov;
pertWin = Parser.Results.pertWin;
initSpeed = Parser.Results.initSpeed;
includeEarlyInitCohPertTime = Parser.Results.includeEarlyInitCohPertTime;
dcpDynCohFile = Parser.Results.dcpDynCohFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
sequences = Parser.Results.sequences;
NumClusters = Parser.Results.NumClusters;
checkUnitType = Parser.Results.checkUnitType;
rateCutoff = Parser.Results.rateCutoff;
winT = Parser.Results.winT;
plotOpts = Parser.Results.plotOpts;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Preallocate
passCutoff = nan(1000,1);
Rinit = nan(1701,3,3,1000);
Rdyn = nan(1701,5,1000);
locations = nan(1000,3);
cellID = nan(1000,100,4);

% ResInit = nan(1701,1000,9*100);
% Minit = nan(1701,1000,3,3,1);
% Ninit = nan(1,1000,3,3,1);
if calcCov.On
    Cinit = nan(length(calcCov.binT),length(calcCov.binT),1000);
else
    Cinit = nan(1,1,1000);
end
VarCEinit = nan(1701,1000);

Ceye = nan([length(winT(1):winT(2)),length(winT(1):winT(2)),length(subjects)]);

indx = 1;

%% Cycle through subjects
for subjecti = 1:length(subjects)
    
    %% Set up subject data
    if ~isempty(dcpObjectFile)
        temp = load(dcpObjectFile{subjecti});
    else
        switch subjects{subjecti}
            case 'ar'
                temp = load('~/Projects/DynamicCoherencePhysiology/ar/dcpObjects/dcpObjects20210406.mat');
            case 'fr'
                temp = load('~/Projects/DynamicCoherencePhysiology/fr/dcpObjects/dcpObjects20230322.mat');
            otherwise
                error('Subject not recognized!')
        end
    end
    dcp = temp.dcp;
    if ~strcmp(dcp{1}.sname,subjects{subjecti})
        error(['dcpObjectFile dcp objects are for ' dcp{1}.sname ', not ' subjects{subjecti}])
    end
    
    if ~isempty(sourceDirectory)
        sourceDirectory = sourceDirectory{subjecti};
    else
        switch subjects{subjecti}
            case 'ar'
                sourceDirectoryTemp = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle';
            case 'fr'
                sourceDirectoryTemp = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Frederick';
            otherwise
                error('Subject not recognized!')
        end
    end
    
    %% Get mean and covariance of each unit from each object file
    [RinitTemp, RdynTemp, cellIDTemp, passCutoffTemp, locationsTemp, ResInitTemp, MinitTemp, NinitTemp, CinitTemp, VarCEinitTemp] = collateFiringRates(dcp,...
        'sourceDirectory',sourceDirectoryTemp,'directions',directions,'chanMap',chanMap{subjecti},...
        'rateCutoff',rateCutoff,'checkUnitType',checkUnitType,...
        'initCohCollate',initCohCollate,'dynCohCollate',dynCohCollate,...
        'calcRes',false,'calcCov',calcCov,'calcVarCE',true);
    
    % Add to data matrices
    Rinit(:,:,:,indx:indx+size(RinitTemp,4)-1) = RinitTemp;
    Rdyn(:,:,indx:indx+size(RdynTemp,3)-1) = RdynTemp;
    cellID(indx:indx+size(cellIDTemp,1)-1,:,1:3) = cellIDTemp;
    cellID(indx:indx+size(cellIDTemp,1)-1,:,4) = subjecti;
    passCutoff(indx:indx+size(passCutoffTemp,1)-1,:) = passCutoffTemp;
    locations(indx:indx+size(locationsTemp,1)-1,:) = locationsTemp;
%     ResInit(:,indx:indx+size(ResInitTemp,2)-1,:) = ResInitTemp;
%     Minit(:,indx:indx+size(MinitTemp,2)-1,:,:,:) = MinitTemp;
%     Ninit(1,indx:indx+size(NinitTemp,2)-1,:,:,:) = NinitTemp;
    Cinit(:,:,indx:indx+size(CinitTemp,3)-1) = CinitTemp;
    VarCEinit(:,indx:indx+size(VarCEinitTemp,2)-1) = VarCEinitTemp;
    
    indx = indx+size(locationsTemp,1);
    
    %% Get residual eye speeds
    [Ceye(:,:,subjecti),~,~] = collateCovariance(dcp,'sourceDirectory',sourceDirectoryTemp,...
        'initCohCollate',initCohCollate,'dynCohCollate',dynCohCollate);
    
end


%%
if initCohCollate
    temp = load([sourceDirectoryTemp '/' dcp{1}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{1}.datapath(end-8:end)]);
    initCoh = temp.initCoh;
end

if dynCohCollate
    filei = 1;
    successfulDynCohLoad = false;
    while ~successfulDynCohLoad
        tempFile = [sourceDirectoryTemp '/' dcpDynCoh.dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
            dcpDynCoh.dcp{filei}.datapath(end-8:end) '.mat'];
        if exist(tempFile,'file')
            temp = load(tempFile);
            successfulDynCohLoad = true;
        else
            filei = filei+1;
        end
    end
    dynCoh = temp.dynCoh;
end


%% 
Rinit = Rinit(:,:,:,1:indx-1);
Rdyn = Rdyn(:,:,1:indx-1);
locations = locations(1:indx-1,:);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);
Cinit = Cinit(:,:,1:indx-1);
VarCEinit = VarCEinit(:,1:indx-1);

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff);
Rdyn = Rdyn(:,:,passCutoff);
locations = locations(passCutoff,:);
cellID = cellID(passCutoff,:,:);
Cinit = Cinit(:,:,passCutoff);
VarCEinit = VarCEinit(:,passCutoff);

%% Remove outlier rates
if initCohCollate
    m = squeeze(max(Rinit,[],[1,2,3]))*1000;    
else
    m = zeros(size(squeeze(max(Rinit,[],[1,2,3]))*1000));
end
if dynCohCollate
    m2 = squeeze(max(Rdyn,[],[1,2]))*1000;
else
    m2 = zeros(size(squeeze(max(Rdyn,[],[1,2]))*1000));
end
Rinit = Rinit(:,:,:,m<=150 & m2<=150);
Rdyn = Rdyn(:,:,m<=150 & m2<=150);
locations = locations(m<=150 & m2<=150,:);
cellID = cellID(m<=150 & m2<=150,:,:);
Cinit = Cinit(:,:,m<=150 & m2<=150);
VarCEinit = VarCEinit(:,m<=150 & m2<=150);

%% Get topology data
functionalTopology = load(FunctionalTopologyFile,'NumClusters','Y','idx','centInd','cellID','subjects');

%% Find VarCE and COV by averaging across residuals for each neuron in a cluster
% for clusti = 1:functionalTopology.NumClusters
%     M(:,functionalTopology.idx==clusti,:,:,:) = sum(...
%         repmat(Ninit(:,functionalTopology.idx==clusti,:,:,:),[size(Minit,1),1,1,1,1]).*...
%         Minit(:,functionalTopology.idx==clusti,:,:,:)/sum(Ninit(1,functionalTopology.idx==clusti,:,:,:),'all'),...
%         [3,4,5]);
%     V = var(ResInit,[],3);
%     FF = V./M;
%     phi = min(FF,[],1);
%     gradVarCE = V - repmat(phi,[size(M,1),1]).*M
% end

%% Average VarCE and Cov matrix by functional
[a,b,c] = intersect(initCoh.neuron_t,calcCov.binT);
for clusti = 1:functionalTopology.NumClusters
    varCE(:,clusti) = mean(VarCEinit(:,functionalTopology.idx==clusti),2);
    Cnew(:,:,clusti) = mean(Cinit(:,:,functionalTopology.idx==clusti),3) - ...
        diag(diag(mean(Cinit(:,:,functionalTopology.idx==clusti),3))) + ...
        diag(mean(VarCEinit(b,functionalTopology.idx==clusti),2));
end

%% Plotting
figure('Name','VarCE by cluster')
plot(initCoh.neuron_t,varCE)
xlabel('Time from motion onset (ms)')
ylabel('VarCE (spikes/s)^2')

figure('Name','CovCE by cluster')
for clusti = 1:functionalTopology.NumClusters
    subplot(3,3,clusti)
    imagesc(calcCov.binT(2:end),calcCov.binT(2:end),Cnew(:,:,clusti))
    xlabel('Time from motion onset (ms)')
    ylabel('Time from motion onset (ms)')
end

