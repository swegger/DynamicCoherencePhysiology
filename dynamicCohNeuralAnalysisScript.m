%% Defaults
binT = 0:100:1300;
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
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/dynCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    if ~isempty(dynCoh.R) || ~any(isnan(dynCoh.R(:))) && size(dynCoh.R,2) > 1
        
        passCutoff = [passCutoff; dynCoh.passCutoff];
        Cnew = nan(length(binT),length(binT),length(dynCoh.unitIndex),length(dir));
        for di = 1:length(dir)
            Cnew(:,:,:,di) = findCovariance(dynCoh,binT,speeds,dir(di),...
                unique(dynCoh.sequences),unique(dynCoh.perturbations));
        end
        
        
        % Select preferred direction
        Rtemp = [];
        Ctemp = [];
        for uniti = 1:length(dynCoh.preferredDirectionRelative)
            ind = find(dir == dynCoh.preferredDirectionRelative(uniti));
            Rtemp = cat(3,Rtemp,dynCoh.R(:,:,uniti,ind));
            Ctemp = cat(3,Ctemp,Cnew(:,:,uniti,ind));
        end
        R = cat(3,R,Rtemp);
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
R = R(:,:,~~passCutoff);
C = C(:,:,~~passCutoff);

%% Remove tail
R = R(dynCoh.neuron_t<=1350,:,:);

%% Remove outlier rates
m = squeeze(max(R,[],[1,2]))*1000;
R = R(:,:,m<=150);
C = C(:,:,m<=150);

%% Mean center
R2 = R;
R2 = R2 - mean(R2,[1,2]);

%% SVD
sz = size(R2);
R3 = reshape(R2,[prod(sz(1:2)),prod(sz(end))]);
[U,S,V] = svd(R3(:,:));

%% PCA
CR = cov(R3);
[vecs,vals] = eig(CR);
vecs = fliplr(vecs);
vals = flipud(diag(vals));
reconstructionN = 5;

for seqi = 1:size(R,2)
    RlowD(:,seqi,:) = permute(vecs(:,1:reconstructionN)'*squeeze(R(:,seqi,:))',[2,3,4,1]);
end

%%
colors = colormap('lines');

%% PCA plotting
figure
controlseq = 5;
for ri = 1:reconstructionN
    subplot(2,reconstructionN,ri)
    for seqi = 1:size(RlowD,2)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),RlowD(:,seqi,ri),'Color',colors(seqi,:))
        hold on
    end
    plotVertical([150+0:300:1500]);
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
    
    subplot(2,reconstructionN,ri+reconstructionN)
    for seqi = 1:size(RlowD,2)
        plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),RlowD(:,seqi,ri)-RlowD(:,controlseq,ri),'Color',colors(seqi,:))
        hold on
    end
    plotVertical([150+0:300:1500]);
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
end

figure
pcdisp = [1 2 3];
subplot(1,2,1)
plot3(RlowD(:,:,pcdisp(1)),RlowD(:,:,pcdisp(2)),RlowD(:,:,pcdisp(3)),'o')
xlabel(['PC ' num2str(pcdisp(1))])
ylabel(['PC ' num2str(pcdisp(2))])
zlabel(['PC ' num2str(pcdisp(3))])
grid on

subplot(1,2,2)
plot3(RlowD(:,:,pcdisp(1))-RlowD(:,5,pcdisp(1)),...
    RlowD(:,:,pcdisp(2))-RlowD(:,5,pcdisp(2)),...
    RlowD(:,:,pcdisp(3))-RlowD(:,5,pcdisp(3)),'o')
xlabel(['PC ' num2str(pcdisp(1))])
ylabel(['PC ' num2str(pcdisp(2))])
zlabel(['PC ' num2str(pcdisp(3))])
grid on

%% Quick estimte of MT input from slip and comparison with modes 3-5
eh = vertcat(dynCoh.eye(:).hvel);
ev = vertcat(dynCoh.eye(:).vvel);
eSpeed = sqrt(eh.^2+ev.^2);
slip = 10 -eSpeed;
figure;
[ax,h1,h2] = plotyy(dynCoh.eye_t,log2(nanmean(slip,1)),...
    dynCoh.neuron_t(dynCoh.neuron_t<=1350),1000*(-RlowD(:,:,3)*vals(3)-RlowD(:,:,4)*vals(4)+RlowD(:,:,5)*vals(5)));
ax(1).YColor = [0 0 0];
h1.Color = [0 0 0];
for seqi = 1:length(seqs)
    h2(seqi).Color = colors(seqi,:);
end
ax(2).YColor = colors(5,:);
axes(ax(1))
xlabel('Time from motion onset (ms)')
ylabel('log_2(slip)')
axes(ax(2))
ylabel('Sum of PCs 3-5')

%% SVD plotting
figure
conds = 1:5;
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

%%
figure
subplot(1,2,1)
mC = mean(C,3);
mCbeh = nanmean(Cbeh,3);
imagesc(binT,binT,mC-diag(diag(mean(C,3))))
axis square
xlabel('Time from motion onset (ms)')
ylabel('Time from motion onset (ms)')

subplot(1,2,2)
% binTbeh = binT(1):binT(end);
% for bi = 1:length(binTbeh)
%     plot(binTbeh(bi+1:end)-binTbeh(bi),mCbeh(bi,bi+1:end)/mCbeh(bi,bi),'Color',[0.8 0.8 0.8])
%     hold on
% end
% for bi = 1:length(binT)
%     plot(binT(bi+1:end)-binT(bi),mC(bi,bi+1:end))
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

xlabel('Time from motion onset (ms)')
ylabel('Normalized autocorrelation')
legend({'Average across neurons','Eye speed'})