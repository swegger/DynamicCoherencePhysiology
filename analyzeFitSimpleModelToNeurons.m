function analyzeFitSimpleModelToNeurons(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs

Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Find files and get data, fit model parameters, ridge regression value

%% Find model result files
switch subject
    case 'ar'
        files = dir('/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle/fitSimpleFEFModelResults/*20231227*.mat');
end

for filei = 1:length(files)
    temp = load([files(filei).folder '/' files(filei).name]);
    if length(temp.modelFEF.tau) == 1
        temp.modelFEF.tau = temp.modelFEF.tau*ones(size(temp.modelFEF.baseLine));
    end
    Rtemp = squeeze(temp.RtoFit(:,:,2,:));
    for neuroni = 1:size(Rtemp,3)
        tempHat = SimpleFEFmodel(temp.modelFEF.W(:,neuroni)',...
            temp.modelFEF.baseLine(neuroni),temp.modelFEF.R0(neuroni),temp.modelFEF.tau(neuroni)/temp.modelFEF.dt,temp.inputs);
        RhatTemp(:,:,neuroni) = tempHat(:,:,2);
    end
    sse(filei) = sum((Rtemp(:) - RhatTemp(:)).^2);
    sseNeuron(:,filei) = sum( (Rtemp - RhatTemp).^2 ,[1,2]);
    lambdaRidge(filei) = temp.lambdaRidge;
    
end

%% Get data from best fitting model over ridge values
[~,minInd] = min(sse);
temp = load([files(minInd).folder '/' files(minInd).name]);
R = temp.RtoFit;
modelFEF = temp.modelFEF;
tau = temp.tau;
dt = temp.dt;
inputs = temp.inputs;
for neuroni = 1:size(R,4)
    Rhat(:,:,:,neuroni) = SimpleFEFmodel(modelFEF.W(:,neuroni)',...
        modelFEF.baseLine(neuroni),modelFEF.R0(neuroni),tau/dt,inputs);
end
t = temp.mt{1}.neuron_t(temp.iMT);
spref = temp.spref;
spref = spref(~isnan(temp.interpolatedR(1,1,1,:)));
[~,prefSort] = sort(spref);

%%

for si = 1:length(speedsFEF)
    rvalAll = corrcoef([squeeze(R(:,si,2,:)) squeeze(Rhat(:,si,2,:))]);
    rvals(:,si) = diag(rvalAll(size(R,4)+1:end,1:size(R,4)));
end

%% Analysis
inputSigma = squeeze(std(temp.interpolatedR,[],[1,2,3]));
inputSigma = inputSigma(~isnan(temp.interpolatedR(1,1,1,:)));

W = modelFEF.W;
Wz = W.*repmat(inputSigma,[1,size(W,2)]);

% Normalize by maximium value of W for each FEF neuron
wmax = max(abs(Wz),[],2);
Weff = Wz./repmat(wmax,[1,size(W,2)]);

% SVD
[U,S,V] = svd(Weff);

% Approach from the important dimensions in Rhat
Y = nan(size(R,1)*size(R,2)*size(R,3),size(R,4));
Yhat = nan(size(Rhat,1)*size(Rhat,2)*size(Rhat,3),size(Rhat,4));
idx = 1;
for si = 1:size(Rhat,2)
    for ci = 1:size(Rhat,3)
        Y(idx:idx+size(R,1)-1,:) = R(:,si,ci,:);
        Yhat(idx:idx+size(Rhat,1)-1,:) = Rhat(:,si,ci,:);
        idx = idx+size(Rhat,1);
    end
end
[Uy,Sy,Vy] = svd(Y,0);
Wy = Wz*Vy;
[Uhat,Shat,Vhat] = svd(Yhat,0);
What = Wz*Vhat;                      % Now ordering weight matrix by the important dimensions in Rhat (e.g. reduced rank regression approach)


%% Plot results

%% Best ridge regression for predicting held-out data
figure
[~,sortInd] = sort(lambdaRidge);
subplot(1,2,1)
loglog(lambdaRidge(sortInd),sse(sortInd),'ko-')
hold on
plotVertical(lambdaRidge(sse == min(sse)));
xlabel('Ridge regression coefficient')
ylabel('Sum of squared errors (spikes/sec)^2')

subplot(1,2,2)
loglog(lambdaRidge(sortInd),sseNeuron(:,sortInd),'-','Color',[0 0 0 0.1])
[~,sseMins] = min(sseNeuron,[],2);
hold on
plotVertical(lambdaRidge(sseMins));
xlabel('Ridge regression coefficient')
ylabel('Sum of squared errors (spikes/sec)^2')


%% Best and worst fits (on average across held-out data)
mrvals = mean(rvals,2);
[~,mrsort] = sort(mrvals);

figure
for ni = 1:5
    subplot(2,5,ni)
    plot(t,R(:,:,2,mrsort(ni))*1000)
    hold on
    plot(t,Rhat(:,:,2,mrsort(ni))*1000,'k')
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    title(['Rank = ' num2str(ni)])
end

for ni = 1:5
    subplot(2,5,ni+5)
    plot(t,R(:,:,2,mrsort(end-ni+1))*1000)
    hold on
    plot(t,Rhat(:,:,2,mrsort(end-ni+1))*1000,'k')
    xlabel('Time from motion onset (ms)')
    ylabel('Spikes/s')
    title(['Rank = ' num2str(length(mrsort)-ni+1)])
end

%% Prediction of held-out data
figure
subplot(1,3,1)
title('Example neuron')
cellID2 = squeeze(temp.cellID(:,1,1:2));
neuroni = find(ismember(cellID2, [77,140],'rows'));
for si = 1:length(speedsFEF)
    plot(t,R(:,si,2,neuroni)*1000,'Color',speedColors(si,:))
    hold on
    plot(t,Rhat(:,si,2,neuroni)*1000,'--','Color',speedColors(si,:))
end
xlabel('Time from motion onset (ms)')
ylabel('Spikes/s')

subplot(1,3,2)
title('Mean across neurons')
for si = 1:length(speedsFEF)
    plot(t,nanmean(R(:,si,2,:),4)*1000,'Color',speedColors(si,:))
    hold on
    plot(t,nanmean(Rhat(:,si,2,:),4)*1000,'--','Color',speedColors(si,:))
end
xlabel('Time from motion onset (ms)')
ylabel('Mean spikes/s (across neurons)')

subplot(1,3,3)
title('R^2 by neuron')
for si = 1:length(speedsFEF)
    plot(rvals(:,si).^2,'o-','Color',speedColors(si,:))
    hold on
end
xlim([1 256])
xlabel('Neuron number')
ylabel('R^2')

%% Model parameters
figure
ah = subplot(1,2,1);
imagesc(modelFEF.W(prefSort,:)')
ah.YDir = 'normal';
xlabel('MT neuron')
ylabel('FEF neuron')

subplot(1,2,2)
plot(modelFEF.baseLine*1000,1:length(modelFEF.baseLine),'o')
axis tight
ylabel('FEF  neuron')
xlabel('Baseline rate')

%% Analysis of weight matrix
figure
subplot(2,2,1)
imagesc(Weff(prefSort,:)')
xlabel('MT neuron')
ylabel('FEF neuron')

subplot(2,2,2)
plot(cumsum(diag(S.^2))./sum(diag(S.^2)),'o')
axis tight
xlabel('Dimensions')
ylabel('Cumulative variance explained')

subplot(2,2,3)
rankN = 1;
for ri = 1:rankN
    semilogx(spref(prefSort),U(prefSort,ri),'o')
    hold on
end
xlim([1 300])
xlabel('Speed preference (deg/s)')
ylabel('MT unit weight')

subplot(2,2,4)
[~,rvalsort] = sort(rvals(:,3));
for ri = 1:rankN
    plot(V(rvalsort,ri),'o')
    hold on
end
xlabel('FEF neuron #')
ylabel('FEF unit weight')

%% Analysis of weight matrix in space that best describes predicted responses
figure
subplot(1,3,1)
plot(cumsum(diag(Shat.^2))./sum(diag(Shat.^2)),'o')
axis tight
xlabel('Dimensions')
ylabel('Cumulative variance explained')

subplot(1,3,2)
rankN = 1;
for ri = 1:rankN
    semilogx(spref(prefSort),What(prefSort,ri),'o')
    hold on
end
xlim([1 300])
xlabel('Speed preference (deg/s)')
ylabel('MT unit weight')

subplot(1,3,3)
[~,rvalsort] = sort(rvals(:,3));
for ri = 1:rankN
    plot(Vhat(rvalsort,ri),'o')
    hold on
end
xlabel('FEF neuron #')
ylabel('FEF unit weight')

figure
ci = 2;
for ri = 1:rankN
    subplot(rankN,4,1+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,-(squeeze(R(:,si,ci,:))*Vhat(:,ri) - modelFEF.baseLine*Vhat(:,ri)),'Color',speedColors(si,:))
        hold on
        plot(t,-(squeeze(Rhat(:,si,ci,:))*Vhat(:,ri) - modelFEF.baseLine*Vhat(:,ri)),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Decay dynamics along dim ' num2str(ri)])
    
    subplot(rankN,4,2+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,What(:,ri)'*inputs(:,:,si,ci),'Color',speedColors(si,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['MT inputs along dim ' num2str(ri)])
    
    subplot(rankN,4,3+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,-(squeeze(R(:,si,ci,:))*Vhat(:,ri) - modelFEF.baseLine*Vhat(:,ri))' + What(:,ri)'*inputs(:,:,si,ci),'Color',speedColors(si,:))
        hold on
        plot(t,-(squeeze(Rhat(:,si,ci,:))*Vhat(:,ri) - modelFEF.baseLine*Vhat(:,ri))' + What(:,ri)'*inputs(:,:,si,ci),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['dr/dt along dim ' num2str(ri)])
    
    subplot(rankN,4,4+(ri-1)*4)
    for si = 1:length(speedsFEF)
        plot(t,squeeze(R(:,si,ci,:))*Vhat(:,ri),'Color',speedColors(si,:))
        hold on
        plot(t,squeeze(Rhat(:,si,ci,:))*Vhat(:,ri),'--','Color',speedColors(si,:))
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['Response along dim ' num2str(ri)])
    
end

%% Functions

%% Simple model
function R = SimpleFEFmodel(W,baseLine,R0,tau,inputs)
    sz = size(inputs);
    R = nan(sz(2:4));
    R(1,:,:) = repmat(R0,[1,sz(3),sz(4)]);
    for ti = 2:sz(2)
        for si = 1:sz(3)
            for ci = 1:sz(4)
                dR(ti,si,ci) = -(R(ti-1,si,ci)-baseLine) + W*inputs(:,ti,si,ci);
                R(ti,si,ci) = R(ti-1,si,ci) + dR(ti,si,ci)/tau;
            end
        end
    end