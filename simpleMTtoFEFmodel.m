function [modelRval,mtWeights,gainPrediction,spref] = simpleMTtoFEFmodel(varargin)
%% mtNeuralAnalysis
%
%
%
%%

%% Defaults

%% Parse input
Parser = inputParser;

addParameter(Parser,'fefN',1000)
addParameter(Parser,'leakRate',NaN)
addParameter(Parser,'fefBaseline',10)
addParameter(Parser,'randWeight',10)
addParameter(Parser,'structuredWeight',0)
addParameter(Parser,'fractTuned',0)
addParameter(Parser,'binT',50:100:850)
addParameter(Parser,'dir',0)
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'pollTimes',[150, 750])
addParameter(Parser,'behaviorFile','neuronTyping20221229.mat')
addParameter(Parser,'tIndex',[2,3,4])
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'weightModel','lowRankDataDriven')
addParameter(Parser,'pltFlg',true)

parse(Parser,varargin{:})

fefN = Parser.Results.fefN;
leakRate = Parser.Results.leakRate;
fefBaseline = Parser.Results.fefBaseline;
randWeight = Parser.Results.randWeight;
structuredWeight = Parser.Results.structuredWeight;
fractTuned = Parser.Results.fractTuned;
binT = Parser.Results.binT;
dir = Parser.Results.dir;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
objectFile = Parser.Results.objectFile;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
pollTimes = Parser.Results.pollTimes;
behaviorFile = Parser.Results.behaviorFile;
tIndex = Parser.Results.tIndex;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
weightModel = Parser.Results.weightModel;
pltFlg = Parser.Results.pltFlg;

if any(isnan(leakRate(:)))
    leakRate = rand(fefN,1);
end

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Load mtObjs
load(objectFile)

%% Get mean of each unit
R = [];
spref = [];
for filei = 1:length(mt)
    disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
    
    Rtemp = nan(length(mt{filei}.neuron_t),length(speeds),length(cohs));
    
    
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            [~,condLogical] = trialSort(mt{filei},0,speeds(si),NaN,cohs(ci));
            Rtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
        end
    end
    
    
    if sprefFromFit
        [mu,sig,~,normR,~] = fitSpeedTuning(mt{filei});
        spref = [spref mu];
        
        if checkMTFit
            s = linspace(min(speeds)-0.1*min(speeds),max(speeds)*1.1,20);
            h = figure;
            semilogx(mt{filei}.speeds,normR,'o')
            hold on
            semilogx(s,speedTuning(mt{filei},s,mu,sig))
            ax = axis;
            text(ax(1),0.95*ax(4),['\mu = ' num2str(mu) ', \sig = ' num2str(sig)])
            input('Press enter to continue ')
            close(h);
        end
    else
        spref = [spref speeds(nansum(Rtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
            max(nansum(Rtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100))))];
    end
    R = cat(4,R,Rtemp);
    %end
end

mtNeuron_t = mt{filei}.neuron_t;


%% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, Ss, Cs, Ss2, Cs2] = interpolateMT_BehlingToEgger(R,speeds,cohs,speedsFEF,cohsFEF);
% [Ss,Cs] = meshgrid(speeds,cohs);
% Ss = Ss';
% Cs = Cs';
% [Ss2,Cs2] = meshgrid(speedsFEF,cohsFEF);
% Ss2 = Ss2';
% Cs2 = Cs2';
% for neuroni = 1:size(R,4)
%     for ti = 1:size(R,1)
%         temp = squeeze(R(ti,:,:,neuroni));
%         if any(sum(~isnan(temp'),1)>2) && any(sum(~isnan(temp'),2)>2)       % Check if there at least two data points on each dimension for interpolation
%             % Strip any row or column that doesnt' have at least 2 data points for interpolation
%             interpolatedR(ti,:,:,neuroni) = interp2(Ss(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
%                 Cs(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
%                 temp(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',Ss2',Cs2','spline')';
%         else
%              interpolatedR(ti,:,:,neuroni) = nan(size(Ss2));
%         end
%     end
% end

%% Load behavioral gain data
beh = load(behaviorFile,'gain','initGain','dynCohT','dynCoh');
for ti = 1:length(pollTimes)
    beh.cohDyn(:,ti) = beh.dynCohT(beh.dynCoh.neuron_t == pollTimes(ti),:);
end

%% Implement a pool of neurons with integration behavior defined by leakRate
dfef = zeros(fefN,size(interpolatedR,1),size(interpolatedR,2),size(interpolatedR,3)); 
fef = zeros(fefN,size(interpolatedR,1),size(interpolatedR,2),size(interpolatedR,3));

switch weightModel
    case 'randomStructuredMix'
        structuredWeight = structuredWeight/mean(log2(spref));
        fef(:,mt{1}.neuron_t<=0,:,:) = fefBaseline;
        mtWeights = randn(fefN,size(R,4));
        mtWeights = (...
            mtWeights*randWeight + ...
            structuredWeight*repmat(log2(spref),[fefN,1]) .* repmat(randsample([-1 1],fefN,true)',[1,size(R,4)])...
            )/fefN;
        speedTunedIndices = randsample(fefN,round(fractTuned*fefN));
        
        mtWeights(speedTunedIndices,:) = randWeight/10*repmat(log2(spref),[length(speedTunedIndices),1]) .* repmat(randsample([-1 1],length(speedTunedIndices),true)',[1,size(R,4)])/fefN;

    case 'lowRankDataDriven'
        rankN = 100;
        for ri = 1:rankN
            rightHand(:,ri) = zeros(size(R,4),1);
            rightHand(spref > 1,ri) = sqrt(log(spref(spref>1).^2))'.*randn(sum(spref>1),1);
            leftHand(:,ri) = randn(fefN,1);
            weightsTemp(:,:,ri) = leftHand(:,ri)*rightHand(:,ri)';
        end
        mtWeights = sum(weightsTemp,3)/rankN/fefN;
               
end

fefTau = 40/(mtNeuron_t(2)-mtNeuron_t(1));
for ti = find(mt{1}.neuron_t==0):length(mt{1}.neuron_t)
    for si = 1:length(speedsFEF)
        for ci = 1:length(cohsFEF)
            mtInputTemp = mtWeights(:,~isnan(interpolatedR(ti,si,ci,:)))*permute(interpolatedR(ti,si,ci,~isnan(interpolatedR(ti,si,ci,:))),[4,1,2,3]);
            dfef(:,ti,si,ci) = -leakRate.*(fef(:,ti-1,si,ci)-fefBaseline) + mtInputTemp;
            fef(:,ti,si,ci) = fef(:,ti-1,si,ci) + dfef(:,ti,si,ci)/fefTau;
            fef(fef < 0) = 0;
        end
    end
end


%% Model analysis
% Regress population against gain
taccept = mtNeuron_t >= 150 & mtNeuron_t <= mtNeuron_t(end);
tempGain = reshape(beh.initGain(:,:,tIndex(2)),[numel(Ss2),1]);
betaFEF = nan(fefN,4);
for neuroni = 1:fefN
    betaFEF(neuroni,:) = regress(reshape(fef(neuroni,taccept,:,:),[numel(Ss2)*sum(taccept),1]),...
        repmat([Ss2(:),Cs2(:),tempGain,ones(numel(Ss2),1)],[sum(taccept),1]));
end
dimLabels = {'Speed','Coherence','Gain','Constant'};

% Find projection of population data along each dimension
[betaFEFOrth,~] = qr(betaFEF);
pBopt = nan(4,size(fef,2),size(fef,3),size(fef,4));
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        pBopt(:,:,si,ci) = betaFEFOrth(:,1:4)'*fef(:,:,si,ci);
    end
end

% Determine how well gain dimension activity predicts behavioral gain
% across time in initCoh task
for di = 1:size(pBopt,1)
    for ti = 1:length(pollTimes)
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                gainPrediction(si,ci,ti+1,di) = pBopt(di,mtNeuron_t == pollTimes(ti),si,ci);
            end
        end
    end
    tempBeh = reshape(beh.initGain(:,:,2:3),[numel(Ss2)*2,1]);
    tempProj = reshape(gainPrediction(:,:,2:3,di),[numel(Ss2)*2,1]);
    rval = corrcoef(tempBeh,tempProj);
    modelRval(di) = rval(1,2);
end

%% Plot

if pltFlg
    
    %% Learky integrator model results
    figure('Name','Model FEFsem activity','Position',[409 375 1914 420])
    samps = randsample(fefN,100);
    for si = 1:length(speedsFEF)
        subplot(1,length(speedsFEF)+1,si)
        fefnorm = (fef(samps,mtNeuron_t>=0,si,2)') ./ ...
            repmat(max(fef(samps,mtNeuron_t>=0,si,2)',[],1),[sum(mtNeuron_t>=0),1]);
        plot(mtNeuron_t(mtNeuron_t>=0),fef(samps,mtNeuron_t>=0,si,2)','Color',[speedColors(si,:) 0.2])
        hold on
        xlabel('Time from motion onset (ms)')
        ylabel('Model FEFsem response')
    end
    subplot(1,length(speedsFEF)+1,length(speedsFEF)+1)
    plot(fef(:,mtNeuron_t == 200,2,2),squeeze(fef(:,mtNeuron_t == 200,1,2)),...
        'o','Color',speedColors(1,:));
    hold on
    plot(fef(:,mtNeuron_t == 200,2,2),squeeze(fef(:,mtNeuron_t == 200,3,2)),...
        'o','Color',speedColors(3,:));
    axis equal
    plotUnity;
    axis square
    xlabel(['Respones to ' num2str(speedsFEF(2)) ' deg/s'])
    ylabel(['Respones to ' num2str(speedsFEF([1,3])) ' deg/s'])

    figure('Name','Projection of regression axes','Position',[405 95 1924 953])
    for si = 1:length(speedsFEF)
        for di = 1:size(pBopt,1)
            subplot(4,length(speedsFEF),si+(di-1)*length(speedsFEF))
            for ci = 1:length(cohsFEF)
                plot(mtNeuron_t(mtNeuron_t>=0),pBopt(di,mtNeuron_t>=0,si,ci),'Color',[speedColors(si,:) cohsFEF(ci)/100])
                hold on
            end
            plotVertical(pollTimes);
            xlabel('Time from motion onset (ms)')
            ylabel(['Projection on ' dimLabels{di} ' axis'])
        end
    end

    figure('Name','Measured gain vs projection on each regression axis','Position',[409 375 1914 420])
    symbs = {'d-','o-'};
    for di = 1:size(pBopt,1)
        subplot(1,size(pBopt,1),di)
        for ti = 1:2
            for si = 1:length(speedsFEF)
                behTemp = beh.initGain(si,:,ti+1);
                projTemp = gainPrediction(si,:,ti+1,di);
                plot(behTemp(:),projTemp(:),symbs{ti},...
                    'Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                hold on
            end
        end
        xlabel('Measured gain')
        ylabel(['Projection on ' dimLabels{di} ' axis'])
        axis square
        ax = axis;
        text(ax(2)*0.9,(ax(4)-ax(3))*0.9+ax(3),['R = ' num2str(modelRval(di))],'HorizontalAlignment','right')
    end
    legend('Early','Late','Location','SouthEast')

    figure('Name','PCA')
    temp = [];
    for si = 1:length(speedsFEF)
        for ci = 1:length(cohsFEF)
            temp = cat(2,temp,fef(:,:,si,ci));
        end
    end
    [COEFF, SCORE] = pca(temp');
    for si = 1:length(speedsFEF)
        for ci = 1:length(cohsFEF)
            A = COEFF(:,:)*fef(:,:,si,ci);
            subplot(1,4,1)
            plot3(A(1,:),A(2,:),A(3,:),'-','Color',[speedColors(si,:) cohsFEF(ci)/100])
            hold on

            for dimi = 1:3
                subplot(1,4,1+dimi)
                plot(mtNeuron_t(mtNeuron_t>=0),A(dimi,mtNeuron_t>=0),'-','Color',[speedColors(si,:) cohsFEF(ci)/100])
                hold on
            end
        end
    end
    subplot(1,4,1)
    grid on
    xlabel('PC_1')
    ylabel('PC_2')
    zlabel('PC_3')
    
end