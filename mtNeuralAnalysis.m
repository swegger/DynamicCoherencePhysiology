function [WA,interpolatedWA,mtNeuron_t,R,interpolatedR,spref,mt] = mtNeuralAnalysis(varargin)
%% mtNeuralAnalysis
%
%
%
%%

%% Defaults

%% Parse input
Parser = inputParser;

addParameter(Parser,'binT',50:100:850)
addParameter(Parser,'dir',[0])
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'pollTimes',[150, 750])
addParameter(Parser,'behaviorFile','neuronTyping20221229.mat')
addParameter(Parser,'tIndex',[2,3,4])
addParameter(Parser,'pltFlg',true)

parse(Parser,varargin{:})

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
pltFlg = Parser.Results.pltFlg;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Load mtObjs
load(objectFile)

%% Get mean and covariance of each unit
C = [];
R = [];
spref = [];
for filei = 1:length(mt)
    disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
    
    Rtemp = nan(length(mt{filei}.neuron_t),length(speeds),length(cohs));
    
    Cnew = nan(length(binT),length(binT),...
        length(mt{filei}.unitIndex),length(dir));
    Cnew = findCovariance(mt{filei},binT,speeds,cohs,dir(1));
    
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            [~,condLogical] = trialSort(mt{filei},0,speeds(si),NaN,cohs(ci));
            Rtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
        end
    end
    %if ~any(isnan(Rtemp(:)))
        spref = [spref speeds(nansum(Rtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
            max(nansum(Rtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100))))];
        R = cat(4,R,Rtemp);
        C = cat(3,C,Cnew);
    %end
end

mtNeuron_t = mt{filei}.neuron_t;

%% Mean center and reshape
R2 = R;
R2 = R2 - nanmean(R2,[1,2]);
sz = size(R2);
R3 = reshape(R2,[prod(sz(1:3)),prod(sz(end))]);
R4 = R3(:,sum(isnan(R3),1) == 0);
R5 = R(:,:,:,sum(isnan(R3),1) == 0);

%% PCA
CR = cov(R4);
[vecs,vals] = eig(CR);
vecs = fliplr(vecs);
vals = flipud(diag(vals));
reconstructionN = 5;

for si = 1:size(R,2)
    for ci = 1:size(R,3)
        RlowD(:,si,ci,:) = permute(vecs(:,1:reconstructionN)'*squeeze(R5(:,si,ci,:))',[2,3,4,1]);
    end
end
%% Feed foward gain model
for si = 1:length(speeds)
    for ci = 1:length(cohs)
%        WA(:,si,ci) = log2(spref)*permute(R(:,si,ci,:),[4,1,2,3]);
        for ti = 1:size(R,1)
            WA(ti,si,ci) = nanmean(log2(spref').*permute(R(ti,si,ci,:),[4,1,2,3]));
        end
    end
end

% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[Ss,Cs] = meshgrid(speeds,cohs);
Ss = Ss';
Cs = Cs';
[Ss2,Cs2] = meshgrid(speedsFEF,cohsFEF);
Ss2 = Ss2';
Cs2 = Cs2';

for ti = 1:size(WA,1)
    temp = WA(ti,:,:);
    interpolatedWA(ti,:,:) = interp2(Ss',Cs',squeeze(temp)',Ss2',Cs2','spline')';
end

%% Interpolate each neuron's responses
for neuroni = 1:size(R,4)
    for ti = 1:size(R,1)
        temp = R(ti,:,:,neuroni);
        if sum(any(isnan(squeeze(temp)),2)) >= size(Ss,2)-1 || sum(any(isnan(squeeze(temp)),1)) >= size(Ss,1)-1
            interpolatedR(ti,:,:,neuroni) = nan(size(Ss2));
        else
            interpolatedR(ti,:,:,neuroni) = interp2(Ss',Cs',squeeze(temp)',Ss2',Cs2','spline')';
        end
    end
end

%% Implement Shunting inhibition similar to Darlington et. al. 2018
WA2 = WA - repmat(WA(1,:,:),[size(WA,1),1,1]); %repmat(mean(WA(mt{1}.neuron_t <= 30,:,:),[1,2,3]),[size(WA,1),1,1]);
WAI = zeros(size(WA2));
WAIfast = zeros(size(WA2));
WAIslow = zeros(size(WA2));
SWA = zeros(size(WA2));
SWAf = zeros(size(WA2));
SWAs = zeros(size(WA2));
tau = 100;
for ti = find(mt{1}.neuron_t==0):length(mt{1}.neuron_t)
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            SWA(ti,si,ci) = WA2(ti,si,ci)./(1+WAI(ti-1,si,ci)*100);
            SWA(SWA<0) = 0;
            dWAI = -0*WAI(ti-1,si,ci) + SWA(ti,si,ci);
            WAI(ti,si,ci) = WAI(ti-1,si,ci) + dWAI/tau;
            WAI(WAI<0) = 0;
            
            
            SWAf(ti,si,ci) = WA2(ti,si,ci)./(1+WAIfast(ti-1,si,ci)*100);
            SWAf(SWA<0) = 0;
            dWAIf = -0.8*WAIfast(ti-1,si,ci) + SWAf(ti,si,ci);
            WAIfast(ti,si,ci) = WAIfast(ti-1,si,ci) + dWAIf/tau;
            WAIfast(WAIfast<0) = 0;
            
            SWAs(ti,si,ci) = WA2(ti,si,ci)./(1+WAIslow(ti-1,si,ci)*100);
            SWAs(SWAs<0) = 0;
            dWAIs = -0.2*WAIslow(ti-1,si,ci) + SWAs(ti,si,ci);
            WAIslow(ti,si,ci) = WAIslow(ti-1,si,ci) + dWAIs/tau;
            WAIslow(WAIslow<0) = 0;
            
        end
    end
end

WAI = WAI - WAI(mt{filei}.neuron_t==0,:,:);
WAIfast = WAIfast - WAIfast(mt{filei}.neuron_t==0,:,:);
WAIslow = WAIslow - WAIslow(mt{filei}.neuron_t==0,:,:);

WAItotal = WAI + WAIfast - WAIslow;

%% Compare to behavior

% Load behavioral gain data
beh = load(behaviorFile,'gain','initGain','dynCohT','dynCoh');
for ti = 1:length(pollTimes)
    beh.cohDyn(:,ti) = beh.dynCohT(beh.dynCoh.neuron_t == pollTimes(ti),:);
end


for ti = 1:length(pollTimes)
    
    % Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
    temp = WA(mt{filei}.neuron_t == pollTimes(ti),:,:);
    temp2 = WAIfast(mt{filei}.neuron_t == pollTimes(ti),:,:);
    gainEstimate(ti,:,:) = interp2(Ss',Cs',squeeze(temp)',Ss2',Cs2','spline')';
    gainEstimate2(ti,:,:) = interp2(Ss',Cs',squeeze(temp2)',Ss2',Cs2','spline');
    
    beta(ti,:) = regress(reshape(gainEstimate(ti,:,:),[numel(Ss2),1]),...
        [reshape(beh.initGain(:,:,tIndex(ti)),[numel(Ss2),1]) ones(numel(Ss2),1)]);
    beta2(ti,:) = regress(reshape(gainEstimate2(ti,:,:),[numel(Ss2),1]),...
        [reshape(beh.initGain(:,:,tIndex(ti)),[numel(Ss2),1]) ones(numel(Ss2),1)]);
end

%% Implement a pool of neurons with diverse integration behavior
fefN = 4000;
leakRate = rand(fefN,1);
fefBaseline = 10;
randWeight = 40;
structuredWeight = 0;
fractTuned = 0;
dfef = zeros(fefN,size(interpolatedWA,1),size(interpolatedWA,2),size(interpolatedWA,3)); 
fef = zeros(fefN,size(interpolatedWA,1),size(interpolatedWA,2),size(interpolatedWA,3));
fef(:,mt{1}.neuron_t<=0,:,:) = fefBaseline;
mtWeights = randn(fefN,size(R,4));
% mtWeights = repmat(log2(spref),[fefN,1]) .* reshape(randsample([-1 1],fefN*size(R,4),true),[fefN,size(R,4)])/fefN;
mtWeights = (...
    mtWeights*randWeight + ...
    structuredWeight*repmat(log2(spref),[fefN,1]) .* repmat(randsample([-1 1],fefN,true)',[1,size(R,4)])...
    )/fefN;
speedTunedIndices = randsample(fefN,round(fractTuned*fefN));
mtWeights(speedTunedIndices,:) = randWeight/10*repmat(log2(spref),[length(speedTunedIndices),1]) .* repmat(randsample([-1 1],length(speedTunedIndices),true)',[1,size(R,4)])/fefN;
fefTau = 40/(mtNeuron_t(2)-mtNeuron_t(1));
for ti = find(mt{1}.neuron_t==0):length(mt{1}.neuron_t)
    for si = 1:length(speedsFEF)
        for ci = 1:length(cohsFEF)
            mtInputTemp = mtWeights(:,~isnan(interpolatedR(ti,si,ci,:)))*permute(interpolatedR(ti,si,ci,~isnan(interpolatedR(ti,si,ci,:))),[4,1,2,3]);
%             mtInputTemp = interpolatedWAtemp(ti,si,ci)/size(R,4);
            dfef(:,ti,si,ci) = -leakRate.*(fef(:,ti-1,si,ci)-fefBaseline) + mtInputTemp;
            fef(:,ti,si,ci) = fef(:,ti-1,si,ci) + dfef(:,ti,si,ci)/fefTau;
            fef(fef < 0) = 0;
        end
    end
end

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


figure('Name','Model FEFsem activity','Position',[409 375 1914 420])
samps = randsample(fefN,100);
for si = 1:length(speedsFEF)
    subplot(1,length(speedsFEF)+1,si)
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

%% Plot

if pltFlg
    %% PCA plotting
    colors = colormap('lines');
    for ri = 1:reconstructionN
        for si = 1:size(RlowD,2)
            subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
            for ci = 1:size(RlowD,3)
                plot(mt{1}.neuron_t,RlowD(:,si,ci,ri),'Color',colors(ci,:))
                hold on
            end
        end
        
        %     subplot(1,reconstructionN,ri)
        %     for si = 1:size(RlowD,2)
        %         plot(initCoh.neuron_t,RlowD(:,si,3,ri),'Color',colors(si,:))
        %         hold on
        %     end
        xlabel('Time from motion onset (ms)')
        ylabel(['PC ' num2str(ri)])
    end
    
    %% Weigthed average response
    speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
    figure('Name','Feedfoward gain','Position',[404 217 1982 420])
    for si = 1:length(speeds)
        subplot(1,length(speeds),si)
        for ci = 1:length(cohs)
            plot(mtNeuron_t,WA(:,si,ci),'Color',(1-cohs(ci)/100)*ones(1,3))%[speedColors(si,:) (cohs(ci))/100])
            hold on
        end
        yls(si,:) = ylim;
    end
    
    for si = 1:length(speeds)
        subplot(1,length(speeds),si)
        ylim([min(yls(:,1)) max(yls(:,2))])
        xlim([mtNeuron_t(1) mtNeuron_t(end)])
        text(mtNeuron_t(end)*0.98,max(yls(:,2))*0.9,['Target speed = ' num2str(speeds(si))],...
            'HorizontalAlignment','right')
        plotVertical(pollTimes);
        xlabel('Time from motion onset (ms)')
        ylabel('Gain estimate (au)')
    end
    
    figure('Name','Feedforward gain vs behavior','Position',[687 421 490*length(pollTimes) 788])
    for ti = 1:length(pollTimes)
        subplot(2,length(pollTimes),ti)
        surf(Ss2,Cs2,squeeze(gainEstimate(ti,:,:)),'EdgeColor','none')
        hold on
%         if pollTimes <= 450
            scatter3(Ss2(:),Cs2(:),beta(ti,1)*reshape(beh.initGain(:,:,tIndex(ti)),[numel(Ss2),1]) + beta(ti,2))
%         end
        
        xlabel('Target speed (deg/s)')
        ylabel('Target coherence (deg/s)')
        zlabel('Gain estimate (a.u.)')
    end
    
    for ti = 1:length(pollTimes)
        subplot(2,length(pollTimes),length(pollTimes)+ti)
        for ci = 1:length(cohsFEF)
%             if pollTimes(ti) == 150 || pollTimes(ti) == 450
                plot(beh.initGain(:,ci,tIndex(ti)),squeeze(gainEstimate(ti,:,ci)),'o-','Color',1-cohsFEF(ci)*ones(1,3)/100)
%             end
            %         if pollTimes(ti) == 450 || pollTimes(ti) == 750 || pollTimes(ti) == 1050
            %             tempEsts = squeeze(gainEstimate(ti,speedsFEF == 10,ci))*ones(size(beh.gain(beh.cohDyn(:,ti) == cohsFEF(ci),tIndexDyn(ti))));
            %
            %             plot(beh.gain(beh.cohDyn(:,ti) == cohsFEF(ci),tIndexDyn(ti)),tempEsts,'o','Color',1-cohsFEF(ci)*ones(1,3)/100)
            %         end
            hold on
        end
        
        xlabel('Behavioral gain')
        ylabel('Gain estimate')
    end
    
    %%
    figure
    subplot(1,2,1)
    mC = mean(C,3);
    imagesc(binT,binT,mC-diag(diag(mean(C,3))))
    axis square
    xlabel('Time from motion onset (ms)')
    ylabel('Time from motion onset (ms)')
    
    subplot(1,2,2)
    for ki = 1:length(binT)
        aC(ki) = mean(diag(mC,ki-1));
    end
    plot(binT-50,aC/aC(1))
    hold on
    xlabel('Time from motion onset (ms)')
    ylabel('Normalized autocorrelation')
    legend({'Average across neurons','Eye speed'})
    
    %%
    figure('Name','Feedfoward gain','Position',[404 217 1982 420])
    for si = 1:length(speeds)
        subplot(1,length(speeds),si)
        for ci = 1:length(cohs)
            plot(mt{filei}.neuron_t(mt{filei}.neuron_t>=0),WAIfast(mt{filei}.neuron_t>=0,si,ci),'Color',[speedColors(si,:) (cohs(ci))/100])
            hold on
        end
        yls(si,:) = ylim;
    end
    
    for si = 1:length(speeds)
        subplot(1,length(speeds),si)
        ylim([min(yls(:,1)) max(yls(:,2))])
        plotVertical(pollTimes);
        xlabel('Time from motion onset (ms)')
        ylabel('Gain estimate (au)')
    end
    
    figure('Name','Feedforward gain vs behavior','Position',[687 421 490*length(pollTimes) 788])
    for ti = 1:length(pollTimes)
        subplot(2,length(pollTimes),ti)
        surf(Ss2,Cs2,squeeze(gainEstimate2(ti,:,:)),'EdgeColor','none')
        hold on
%         if pollTimes <= 450
            scatter3(Ss2(:),Cs2(:),beta2(ti,1)*reshape(beh.initGain(:,:,tIndex(ti)),[numel(Ss2),1]) + beta(ti,2))
%         end
        
        xlabel('Target speed (deg/s)')
        ylabel('Target coherence (deg/s)')
        zlabel('Gain estimate (a.u.)')
    end
    
    for ti = 1:length(pollTimes)
        subplot(2,length(pollTimes),length(pollTimes)+ti)
        for ci = 1:length(cohsFEF)
%             if pollTimes(ti) == 150 || pollTimes(ti) == 450
                plot(beh.initGain(:,ci,tIndex(ti)),squeeze(gainEstimate2(ti,:,ci)),'o-','Color',1-cohsFEF(ci)*ones(1,3)/100)
%             end
            %         if pollTimes(ti) == 450 || pollTimes(ti) == 750 || pollTimes(ti) == 1050
            %             tempEsts = squeeze(gainEstimate(ti,speedsFEF == 10,ci))*ones(size(beh.gain(beh.cohDyn(:,ti) == cohsFEF(ci),tIndexDyn(ti))));
            %
            %             plot(beh.gain(beh.cohDyn(:,ti) == cohsFEF(ci),tIndexDyn(ti)),tempEsts,'o','Color',1-cohsFEF(ci)*ones(1,3)/100)
            %         end
            hold on
        end
        
        xlabel('Behavioral gain')
        ylabel('Gain estimate')
    end
    
    %% Example MT neurons
    figure('Name','Example MT neuron PSTHs','Position',[329 736 2197 586])
    exInds = [1,5];
    
    ylims = [0 0];
    for exNeuroni = 1:length(exInds)
        for si = 1:length(speeds)
            subplot(length(exInds),length(speeds),si+(exNeuroni-1)*length(speeds))
            for ci = 1:length(cohs)
                plot(mtNeuron_t,R(:,si,ci,exInds(exNeuroni)),...
                    'Color',(1-cohs(ci)/100)*ones(1,3))
                hold on
            end
            temp_lims = ylim;
            if temp_lims(1) < ylims(1)
                ylims(1) = temp_lims(1);
            end
            if temp_lims(2) > ylims(2)
                ylims(2) = temp_lims(2);
            end
        end
    end
    
    for exNeuroni = 1:length(exInds)
        for si = 1:length(speeds)
            subplot(length(exInds),length(speeds),si+(exNeuroni-1)*length(speeds))
            ylim(ylims);
            xlim([mtNeuron_t(1),mtNeuron_t(end)])
            text(mtNeuron_t(end)*0.98,ylims(2)*0.8,['Target speed = ' num2str(speeds(si))],...
                'HorizontalAlignment','right')
            xlabel('Time from motion onset (ms)')
            ylabel('Spikes/s')
        end
    end
               
    
end