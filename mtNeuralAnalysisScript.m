%% Defaults
binT = 50:100:850;
dir = [0];
speeds = [2,4,8,16,32];
cohs = [10 30 70 100];

%% Load mtObjs
load('arMT_20210823.mat')

%% Get mean and covariance of each unit
C = [];
R = [];
for filei = 1:length(mt)
    disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
    
    Cnew = nan(length(binT),length(binT),...
        length(mt{filei}.unitIndex),length(dir));
        Cnew = findCovariance(mt{filei},binT,speeds,cohs,dir(1));
    
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            [~,condLogical] = trialSort(mt{filei},0,speeds(si),NaN,cohs(ci));
            Rtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
        end
    end
    R = cat(4,R,Rtemp);
    C = cat(3,C,Cnew);
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


%%
colors = colormap('lines');

%% PCA plotting
figure
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