function MTtoFEFmodelInstantiations_analysis(subject,varargin)
%%

%%

%% Defaults


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'keepWeights',false)
addParameter(Parser,'behaviorFile','neuronTyping20221229.mat')
addParameter(Parser,'pollTimes',[150, 750])

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
keepWeights = Parser.Results.keepWeights;
behaviorFile = Parser.Results.behaviorFile;
pollTimes = Parser.Results.pollTimes;

%% Find model result files
switch subject
    case 'ar'
        files = dir('/home/seth/Projects/DynamicCoherencePhysiology/ar/MTtoFEFmodelResults/modelRun20231204*.mat');
end

%% Preallocate
temp = load(files(1).name);
fefN = temp.fefN;
rvals = nan(length(files),4);
if keepWeights
    mtWeights = nan(size(temp.mtWeights,1),size(temp.mtWeights,2),length(files));
end
weights.det = nan(length(files),1);
weights.dist = nan(length(files),1);
structuredWeights = repmat(log2(logspace(-1,2,188)),[fefN,1]);

%% Get results from each file
for filei = 1:length(files)
    disp(filei)
    temp = load(files(filei).name);
    rvals(filei,:) = temp.modelRval;
    if keepWeights
        mtWeights(:,:,filei) = temp.mtWeights;
    end
    
    % Find determinant and distance metrics
    weights.det(filei) = det(temp.mtWeights'*temp.mtWeights);
    weights.dist(filei) = norm(reshape(temp.mtWeights,[size(temp.mtWeights,1)*size(temp.mtWeights,2) 1])-structuredWeights(:));
end

%%
contrastR = (rvals(:,3).^2-rvals(:,1).^2)./(rvals(:,3).^2 + rvals(:,1).^2);
[~,minInd] = min(contrastR);
[~,maxInd] = max(contrastR);
temp = load(files(minInd).name);
minW = temp.mtWeights;

temp = load(files(maxInd).name);
maxW = temp.mtWeights;

%% Load behavioral gain data
beh = load(behaviorFile,'gain','initGain','dynCohT','dynCoh');
for ti = 1:length(pollTimes)
    beh.cohDyn(:,ti) = beh.dynCohT(beh.dynCoh.neuron_t == pollTimes(ti),:);
end

%% Plot results

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',size(temp.gainPrediction,1),'sampleN',size(temp.gainPrediction,1));
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%%
figure
plot(rvals(:,1).^2,rvals(:,3).^2,'o')
hold on
axis equal
plotUnity;
axis square

%%
figure
for vali = 1:size(rvals,2)
    subplot(1,size(rvals,2),vali)
    scatter(weights.det*1e20,weights.dist,20,rvals(:,1).^2)
    axis tight
    xlabel('Matrix autocorrelation determinant')
    ylabel('Distance from structured weight matrix')
end

%%
temp = load(files(maxInd).name);
figure('Name','Measured gain vs projection on each regression axis','Position',[409 375 1914 420])
symbs = {'d-','o-'};
for di = 1:size(rvals,2)
    subplot(1,size(rvals,2),di)
    for ti = 1:2
        for si = 1:size(temp.gainPrediction,1)
            behTemp = beh.initGain(si,:,ti+1);
            projTemp = temp.gainPrediction(si,:,ti+1,di);
            plot(behTemp(:),projTemp(:),symbs{ti},...
                'Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
            hold on
        end
    end
    xlabel('Measured gain')
    ylabel(['Projection on ' num2str(di) ' axis'])
    axis square
    ax = axis;
    text(ax(2)*0.9,(ax(4)-ax(3))*0.9+ax(3),['R = ' num2str(rvals(maxInd,di))],'HorizontalAlignment','right')
end
legend('Early','Late','Location','SouthEast')
