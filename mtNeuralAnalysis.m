function mtNeuralAnalysis(varargin)
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
addParameter(Parser,'objectFile','arMT_20210823.mat')
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'pollTimes',[150, 450, 750])
addParameter(Parser,'behaviorFile','neuronTyping20221229.mat')
addParameter(Parser,'tIndex',[2,3,4])

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
    if ~any(isnan(Rtemp(:)))
        spref = [spref speeds(Rtemp(mt{filei}.neuron_t == 50,:,cohs == 100) == max(Rtemp(mt{filei}.neuron_t == 50,:,cohs == 100)))];
        R = cat(4,R,Rtemp);
        C = cat(3,C,Cnew);
    end
end

%% Mean center and reshape
R2 = R;
R2 = R2 - mean(R2,[1,2]);
sz = size(R2);
R3 = reshape(R2,[prod(sz(1:3)),prod(sz(end))]);

%% PCA
CR = cov(R3);
[vecs,vals] = eig(CR);
vecs = fliplr(vecs);
vals = flipud(diag(vals));
reconstructionN = 5;

for si = 1:size(R,2)
    for ci = 1:size(R,3)
        RlowD(:,si,ci,:) = permute(vecs(:,1:reconstructionN)'*squeeze(R(:,si,ci,:))',[2,3,4,1]);
    end
end

%% Feed foward gain model
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        WA(:,si,ci) = log2(spref)*permute(R(:,si,ci,:),[4,1,2,3]);
    end
end

% Load behavioral gain data
beh = load(behaviorFile,'gain','initGain','dynCohT','dynCoh');
for ti = 1:length(pollTimes)
    beh.cohDyn(:,ti) = beh.dynCohT(beh.dynCoh.neuron_t == pollTimes(ti),:);
end

% Interpolate to speed/coh response for dynCoh and initCoh data at
% appropriate time points
[Ss,Cs] = meshgrid(speeds,cohs);
Ss = Ss';
Cs = Cs';
[Ss2,Cs2] = meshgrid(speedsFEF,cohsFEF);
Ss2 = Ss2';
Cs2 = Cs2';
for ti = 1:sum(pollTimes<=450)
    temp = WA(mt{filei}.neuron_t == pollTimes(ti),:,:);
    gainEstimate(ti,:,:) = interp2(Ss',Cs',squeeze(temp)',Ss2',Cs2','spline')';
    
    beta(ti,:) = regress(reshape(gainEstimate(ti,:,:),[numel(Ss2),1]),...
        [reshape(beh.initGain(:,:,tIndex(ti)),[numel(Ss2),1]) ones(numel(Ss2),1)]);
end
            



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
figure('Name','Feedfoward gain','Position',[404 217 1982 420])
for si = 1:length(speeds)
    subplot(1,length(speeds),si)
    for ci = 1:length(cohs)
        plot(mt{filei}.neuron_t,WA(:,si,ci))
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
for ti = 1:sum(pollTimes<=450)
    subplot(2,length(pollTimes),ti)
    surf(Ss2,Cs2,squeeze(gainEstimate(ti,:,:)),'EdgeColor','none')
    hold on
    if pollTimes <= 450
        scatter3(Ss2(:),Cs2(:),beta(ti,1)*reshape(beh.initGain(:,:,tIndex(ti)),[numel(Ss2),1]) + beta(ti,2)) 
    end
    
    xlabel('Target speed (deg/s)')
    ylabel('Target coherence (deg/s)')
    zlabel('Gain estimate (a.u.)')
end

for ti = 1:length(pollTimes)
    subplot(2,length(pollTimes),length(pollTimes)+ti)
    for ci = 1:length(cohsFEF)
        if pollTimes(ti) == 150 || pollTimes(ti) == 450
            plot(beh.initGain(:,ci,tIndex(ti)),squeeze(gainEstimate(ti,:,ci)),'o-','Color',1-cohsFEF(ci)*ones(1,3)/100)
        end
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