%% Defaults
binT = 0:100:1400;
dir = [0 180];

%% Get mean and covariance of each unit
load('dcpObjects20210406.mat');
sourceDirectory = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle';

C = [];
passCutoff = [];
R = [];
Cbeh = [];
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    if ~isempty(initCoh.R)
        
        passCutoff = [passCutoff; initCoh.passCutoff];
        Cnew = nan(length(binT),length(binT),length(initCoh.unitIndex),length(dir));
        for di = 1:length(dir)
            Cnew(:,:,:,di) = findCovariance(initCoh,binT,speeds,cohs,dir(di));
        end
                
        % Select preferred direction
        Rtemp = [];
        Rtemp = [];
        Ctemp = [];
        for uniti = 1:length(initCoh.preferredDirectionRelative)
            ind = find(dir == initCoh.preferredDirectionRelative(uniti));
            Rtemp = cat(4,Rtemp,initCoh.R(:,:,:,uniti,ind));
            Ctemp = cat(3,Ctemp,Cnew(:,:,uniti,ind));
        end
        R = cat(4,R,Rtemp);
        C = cat(3,C,Ctemp);
        
        if strcmp('20210326a',dcp{filei}.datapath(end-8:end))
            disp(['Problems with eye trace for data associated with 20210326a; skipping behavioral covariance...'])
        else
            CbehTemp = findCovarianceBehavior(initCoh,unique(initCoh.speeds),unique(initCoh.coh),0,[0 1400]);
            Cbeh = cat(3,Cbeh,CbehTemp);
        end
    end
end

%% Remove data that doesn't pass cutoff
R = R(:,:,:,~~passCutoff,:);
C = C(:,:,~~passCutoff);

%% Remove outlier rates
m = squeeze(max(R,[],[1,2,3,5]))*1000;
R = R(:,:,:,m<=150,:);
C = C(:,:,m<=150);

%% Mean center
R2 = R(:,:,:,:,1);
R2 = R2 - mean(R2,[1,2,3]);

%% SVD
sz = size(R2);
R3 = reshape(R2,[prod(sz(1:3)),prod(sz(end))]);
[U,S,V] = svd(R3(:,:));

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
    subplot(1,reconstructionN,ri)
%     for ci = 1:size(RlowD,3)
%         plot(initCoh.neuron_t,RlowD(:,2,ci,ri),'Color',colors(ci,:))
%         hold on
%     end
    for si = 1:size(RlowD,2)
        plot(initCoh.neuron_t,RlowD(:,si,3,ri),'Color',colors(si,:))
        hold on
    end
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
end

%% SVD plotting
figure
subplot(1,3,1)
conds = [1 4 7];
ind = 0;
for condi = conds
    ind = ind+1;
    plot3(U((condi-1)*sz(1)+1:condi*sz(1),1),...
        U((condi-1)*sz(1)+1:condi*sz(1),2),...
        U((condi-1)*sz(1)+1:condi*sz(1),3),...
        '-o','Color',colors(ind,:))
    hold on
    
    plot3(U((condi-1)*sz(1)+1,1),...
        U((condi-1)*sz(1)+1,2),...
        U((condi-1)*sz(1)+1,3),...
        'o','Color',colors(ind,:),'MarkerFaceColor',[0 0 0],'MarkerSize',10)
end
grid on
xlabel('1')
ylabel('2')
zlabel('3')
view([90 15])

subplot(1,3,2)
ind = 0;
conds = [2 5 8];
for condi = conds
    ind = ind+1;
    plot3(U((condi-1)*sz(1)+1:condi*sz(1),1),...
        U((condi-1)*sz(1)+1:condi*sz(1),2),...
        U((condi-1)*sz(1)+1:condi*sz(1),3),...
        '-o','Color',colors(ind,:))
    hold on
    
    plot3(U((condi-1)*sz(1)+1,1),...
        U((condi-1)*sz(1)+1,2),...
        U((condi-1)*sz(1)+1,3),...
        'o','Color',colors(ind,:),'MarkerFaceColor',[0 0 0],'MarkerSize',10)
end
grid on
xlabel('1')
ylabel('2')
zlabel('3')
view([90 15])

subplot(1,3,3)
conds = [3 6 9];
ind = 0;
for condi = conds
    ind = ind+1;
    plot3(U((condi-1)*sz(1)+1:condi*sz(1),1),...
        U((condi-1)*sz(1)+1:condi*sz(1),2),...
        U((condi-1)*sz(1)+1:condi*sz(1),3),...
        '-o','Color',colors(ind,:))
    hold on
    
    plot3(U((condi-1)*sz(1)+1,1),...
        U((condi-1)*sz(1)+1,2),...
        U((condi-1)*sz(1)+1,3),...
        'o','Color',colors(ind,:),'MarkerFaceColor',[0 0 0],'MarkerSize',10)
end
grid on
xlabel('1')
ylabel('2')
zlabel('3')
view([90 15])

%%
figure
subplot(1,2,1)
mC = mean(C,3);
mCbeh = nanmean(Cbeh,3);
imagesc(binT,binT,mC-diag(diag(mean(C,3))))
axis square

subplot(1,2,2)
% binTbeh = binT(1):binT(end);
% for bi = 1:length(binTbeh)
%     plot(binTbeh(bi+1:end)-binTbeh(bi),mCbeh(bi,bi+1:end)/mCbeh(bi,bi),'Color',[0.8 0.8 0.8])
%     hold on
% end
% for bi = 1:length(binT)
%     plot(binT(bi+1:end)-binT(bi),mC(bi,bi+1:end)/mC(bi,bi),'Color',[0 0 0],'LineWidth',2)
%     hold on
% end
for ki = 1:length(binT)
    aC(ki) = mean(diag(mC,ki-1));
end
plot(binT,aC/aC(1))
hold on
for ki = 1:size(mCbeh,2)
    aCbeh(ki) = mean(diag(mCbeh,ki-1));
end
plot(binTbeh,aCbeh/aCbeh(1))