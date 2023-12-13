function modelFEF = fitSimpleFEFmodelToNeurons(dcp,varargin)
%% neuronPartialCorrelation
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'P0',NaN)
addParameter(Parser,'ub',NaN)
addParameter(Parser,'lb',NaN)
addParameter(Parser,'tWin',[0 900])
addParameter(Parser,'tau',20)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
directions = Parser.Results.directions;
P0 = Parser.Results.P0;
ub = Parser.Results.ub;
lb = Parser.Results.lb;
tWin = Parser.Results.tWin;
tau = Parser.Results.tau;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Get firing rate data
passCutoff = nan(1000,1);
Rinit = nan(1701,3,3,1000);
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
        
        passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;
        
        % Get data for each neuron
        for uniti = 1:length(initCoh.preferredDirectionRelative)
            ind = find(directions == initCoh.preferredDirectionRelative(uniti));
            Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
            
            
            for j = 1:length(initCoh.unitIndex)
                cellID(indx,j,1) = filei;
                cellID(indx,j,2) = initCoh.unitIndex(uniti);
                cellID(indx,j,3) = initCoh.unitIndex(j);
            end
                        
            indx = indx+1;
        end
    end
end


Rinit = Rinit(:,:,:,1:indx-1);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);

taccept = initCoh.neuron_t >= tWin(1) & initCoh.neuron_t <= tWin(2);

Rinit = Rinit(taccept,:,:,:);

fef_t = initCoh.neuron_t(taccept);

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff);
cellID = cellID(passCutoff,:,:);

%% Get MT data and organize as an input

% Load mtObjs
load(objectFile)

% Get mean of each unit
MT = [];
spref = [];
for filei = 1:length(mt)
    disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
    
    MTtemp = nan(length(mt{filei}.neuron_t),length(speeds),length(cohs));
    
    
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            [~,condLogical] = trialSort(mt{filei},0,speeds(si),NaN,cohs(ci));
            MTtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
        end
    end
    %if ~any(isnan(Rtemp(:)))
        spref = [spref speeds(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
            max(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100))))];
        MT = cat(4,MT,MTtemp);
    %end
end

mtNeuron_t = mt{filei}.neuron_t;


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);
% [Ss,Cs] = meshgrid(speeds,cohs);
% Ss = Ss';
% Cs = Cs';
% [Ss2,Cs2] = meshgrid(speedsFEF,cohsFEF);
% Ss2 = Ss2';
% Cs2 = Cs2';
% for neuroni = 1:size(MT,4)
%     for ti = 1:size(MT,1)
%         temp = squeeze(MT(ti,:,:,neuroni));
%         if any(sum(~isnan(temp'),1)>2) && any(sum(~isnan(temp'),2)>2) % Check if there at least two data points on each dimension for interpolation
%             % Strip any row or column that doesnt' have at least 2 data points for interpolation
%             interpolatedR(ti,:,:,neuroni) = interp2(Ss(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
%                 Cs(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
%                 temp(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',Ss2',Cs2','spline')';
%         else
%              interpolatedR(ti,:,:,neuroni) = nan(size(Ss2));
%         end
%     end
% end

inputs = permute(interpolatedR,[4,1,2,3]);
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);

%% Iterate across neurons and fit data

% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,mtNeuron_t);
RtoFit = Rinit(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);


%% Fit dynamic model
if any(isnan(P0))
    P0 = [0.10,randn(1,size(inputs,1))/1000];
end
if any(isnan(lb))
    lb = [0, -Inf*ones(1,size(inputs,1))];
end
if any(isnan(ub))
    ub = [100, Inf*ones(1,size(inputs,1))];
end
OPTIONS = optimset(@fmincon);
OPTIONS.MaxFunEvals = 3e10;
for neuroni = 1:size(RtoFit,4)
    tic
    disp(['Neuron ' num2str(neuroni) ' of ' num2str(size(RtoFit,4))])
    Rtemp = RtoFit(:,:,[1, length(cohsFEF)],neuroni);
    R0 = nanmean(Rtemp(1,:,:),[2,3]);
    minimizer = @(P)minimizant(P,R0,tau,inputs(:,:,:,[1,length(cohsFEF)]),Rtemp);
    [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub,[],OPTIONS);
    
    modelFEF.R0(neuroni) = R0;
    modelFEF.baseLine(neuroni) = P(1);
    modelFEF.W(:,neuroni) = P(2:end);
    modelFEF.fval(:,neuroni) = fval;
    modelFEF.exitflag(:,neuroni) = exitflag;
    modelFEF.output(:,neuroni) = output;
    modelFEF.lambda(:,neuroni) = lambda;
    modelFEF.grad(:,neuroni) = grad;
    modelFEF.hessian(:,:,neuroni) = hessian;
    toc
end

%% Regression models

% Build matrices
Xfit = nan(size(inputs,2)*size(inputs,3)*(size(inputs,4)-1),size(inputs,1));
Yfit = nan(size(RtoFit,1)*size(RtoFit,2)*(size(RtoFit,3)-1),size(RtoFit,4));

idx = [1, 1];
for si = 1:length(speedsFEF)
    for ci = [1 length(cohsFEF)]
        Xfit(idx(1):idx(1)+size(inputs,2)-1,:) = permute(inputs(:,:,si,ci),[2,1]);
        Yfit(idx(2):idx(2)+size(RtoFit,1)-1,:) = squeeze(RtoFit(:,si,ci,:));
        idx(1) = idx(1) + size(inputs,2);
        idx(2) = idx(2) + size(RtoFit,1);
    end
end

Xtest = nan(size(inputs,2)*size(inputs,3),size(inputs,1));
Ytest = nan(size(RtoFit,1)*size(RtoFit,2),size(RtoFit,4));
idx = [1, 1];
for si = 1:length(speedsFEF)
    Xtest(idx(1):idx(1)+size(inputs,2)-1,:) = permute(inputs(:,:,si,2),[2,1]);
    Ytest(idx(2):idx(2)+size(RtoFit,1)-1,:) = squeeze(RtoFit(:,si,2,:));
    idx(1) = idx(1) + size(inputs,2);
    idx(2) = idx(2) + size(RtoFit,1);
end

% Regress
ridgeLambda = logspace(-1,12,10);
trainBeta = nan([size(Xfit,2),size(Yfit,2),length(ridgeLambda)]);
for ri = 1:length(ridgeLambda)
    trainBeta(:,:,ri) = inv(Xfit'*Xfit + ridgeLambda(ri)*eye(size(Xfit'*Xfit)))*Xfit'*Yfit;
    trainYfit = Xfit*trainBeta(:,:,ri);
    testYfit = Xtest*trainBeta(:,:,ri);
    
    fitQ(ri,1) = sum( (Yfit(:) - trainYfit(:)).^2 );
    fitQ(ri,2) = sum( (Ytest(:) - testYfit(:)).^2 );
end
Beta = trainBeta(:,:,fitQ(:,2)==min(fitQ(:,2)));

% Predicted responses
YhatTrain = Xfit*Beta;
YhatTest = Xtest*Beta;

% Reduced rank regression
rankN = 80;
[Uhat,Shat,Vhat] = svd(YhatTrain,0);
BetaRRR = Beta*Vhat(:,1:rankN)*Vhat(:,1:rankN)';
BetaNull = Beta*Vhat(:,end)*Vhat(:,end)';
YrrrTrain = Xfit*BetaRRR;
YnullTrain = Xfit*BetaNull;
YrrrTest = Xtest*BetaRRR;
YnullTest = Xtest*BetaNull;

Rhat = nan(size(RtoFit));
Rnull = nan(size(RtoFit));
idx = 1;
idxTest = 1;
for si = 1:length(speedsFEF)
    for ci = [1, length(cohsFEF)]
        Rhat(:,si,ci,:) = YrrrTrain(idx:idx+size(RtoFit,1)-1,:);
        Rnull(:,si,ci,:) = YnullTrain(idx:idx+size(RtoFit,1)-1,:);
        idx = idx+size(RtoFit,1);
    end
    
    ci = 2;
    Rhat(:,si,ci,:) = YrrrTest(idxTest:idxTest+size(RtoFit,1)-1,:);
    Rnull(:,si,ci,:) = YnullTest(idxTest:idxTest+size(RtoFit,1)-1,:);
    idxTest = idxTest+size(RtoFit,1);
end

% Find correlation with held out data and do Fisher transformation
rvals = corrcoef([Ytest YrrrTest]);
RhatCC = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZ = 0.5*log( (1+rhatvals)./(1-rhatvals) );


%% Analysis

% Sort MT dimension by speed preference, remove poor fit FEF units, and
% remove leak rate
[~,sortInd] = sort(spref(~isnan(interpolatedR(1,1,1,:))));
spref2 = spref(~isnan(interpolatedR(1,1,1,:)));
spref2 = spref2(sortInd);
cutThreshold = quantile(modelFEF.fval,0.975);
W = modelFEF.W(sortInd,modelFEF.fval<cutThreshold)' ./ ...
    repmat(modelFEF.leakRate(modelFEF.fval<cutThreshold)',[1,size(modelFEF.W,1)]);

% Normalize by maximium value of W for each FEF neuron
wmax = max(W,[],2);
Weff = W./repmat(wmax,[1,size(W,2)]);

% SVD
[U,S,V] = svd(Weff);

%% Plotting

%% Fit quality
figure
for neuroni = size(RtoFit,4)
    tempFEF = RtoFit(:,:,:,neuroni);
    test = SimpleFEFmodel(modelFEF.W(:,neuroni)',modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),tau,inputs);
    for ci = 1:length(cohsFEF)
        subplot(1,length(cohsFEF),ci)
        for si = 1:length(speedsFEF)
            plot(tempFEF(:,si,ci)*1000,test(:,si,ci)*1000,'o','Color',speedColors(si,:))
            hold on
        end
    end
end

for ci = 1:length(cohsFEF)
    subplot(1,length(cohsFEF),ci)
%     axis(1000*[min(RtoFit(:))-0.1*min(RtoFit(:)), max(RtoFit(:))+0.01*max(RtoFit(:)), ...
%         min(RtoFit(:))-0.1*min(RtoFit(:)), max(RtoFit(:))+0.01*max(RtoFit(:))])
    axis square
    plotUnity;
    xlabel('FEF data (sp/s)')
    ylabel('Fit data (sp/s)')
end

%% Reduced rank regression analysis
figure
for ci = 1:length(cohsFEF)
    subplot(1,length(cohsFEF),ci)
    for si = 1:length(speedsFEF)
        plot(squeeze(RtoFit(:,si,ci,:))*1000,squeeze(Rhat(:,si,ci,:))*1000,'o','Color',speedColors(si,:))
        hold on
    end
    plotUnity;
    xlabel('FEF data (sp/s)')
    ylabel('RRR prediction (sp/s)')
    
end

figure
plot(RhatZ,'o')
hold on
plotHorizontal(2.34/sqrt(size(Ytest,1)-size(Xfit,2)));       % Approximate pval of 0.01
xlabel('FEF neuron #')
ylabel('Fisher transformed r-values')

%%
figure('Name','MT output channels found by RRR','Position',[674 218 1837 1104])
rankN = 5;
for ri = 1:rankN
    for ci = 1:length(cohsFEF)
        subplot(rankN,length(cohsFEF),ci + (ri-1)*length(cohsFEF))
        for si = 1:length(speedsFEF)
            MToutputChan(:,si,ci,ri) = inputs(:,:,si,ci)'*Beta*Vhat(:,ri);
            plot(MToutputChan(:,si,ci,ri),'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (steps)')
        ylabel(['Output along dimension#' num2str(ri)])
        tempLims(ci,:) = ylim;
    end
    
    for ci = 1:length(cohsFEF)
        subplot(rankN,length(cohsFEF),ci + (ri-1)*length(cohsFEF))
        ylim([min(tempLims(:,1)) max(tempLims(:,2))])
    end
end

%% Functions

%% Function to fit
function R = SimpleFEFmodel(W,baseLine,R0,tau,inputs)
    sz = size(inputs);
    R = nan(sz(2:4));
    R(1,:,:) = repmat(R0,[1,sz(3),sz(4)]);
    for ti = 2:sz(2)
        for si = 1:sz(3)
            for ci = 1:sz(4)
                dR(ti,si,ci) = (R(ti-1,si,ci)-baseLine) + W*inputs(:,ti,si,ci);
                R(ti,si,ci) = R(ti-1,si,ci) + dR(ti,si,ci)/tau;
            end
        end
    end
    
 %% Objective
function out = minimizant(P,R0,tau,inputs,R)
    baseLine = P(1);
    W = P(2:end);
    Rest = SimpleFEFmodel(W,baseLine,R0,tau,inputs);
    out = sum((R(:) - Rest(:)).^2);