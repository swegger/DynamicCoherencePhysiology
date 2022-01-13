function initiateCohNeuralAnalysisScript(varargin)
%%
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'dcpObjectFile','dcpObjects20210406.mat')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'binT',0:100:1400)
addParameter(Parser,'dir',[0 180])
addParameter(Parser,'initWin',[150 200])
addParameter(Parser,'chanMap',[7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7; 0 8; 23 9; 22 10; 21 11; 20 12; 19 13; 18 14; 17 15; 16 16; 15 17; 14 18; 13 19; 12 20; 11 21; 10 22; 9 23; 8 24])
addParameter(Parser,'win',-200:200)
addParameter(Parser,'calcCC',false)
addParameter(Parser,'calcRSC',true)
addParameter(Parser,'shuffleN',100)
addParameter(Parser,'sampleRate',40)
addParameter(Parser,'ClusterMethod','densityClust')

parse(Parser,varargin{:})

dcpObjectFile = Parser.Results.dcpObjectFile;
sourceDirectory = Parser.Results.sourceDirectory;
binT = Parser.Results.binT;
dir = Parser.Results.dir;
initWin = Parser.Results.initWin;
chanMap = Parser.Results.chanMap;
win = Parser.Results.win;
calcCC = Parser.Results.calcCC;
calcRSC = Parser.Results.calcRSC;
shuffleN = Parser.Results.shuffleN;
sampleRate = Parser.Results.sampleRate;
ClusterMethod = Parser.Results.ClusterMethod;

%% Load dcp object file
load(dcpObjectFile);

%% Get mean and covariance of each unit
C = [];
passCutoff = [];
R = [];
Rsplit = [];
Cbeh = [];
locations = [];
cellID = nan(400,100,3);
CC = nan(400,100,length(win));
CC_CI = nan(400,100,2);
RSC = nan(400,100);
RSC_pval = nan(400,100);
indx = 1;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    if ~isempty(initCoh.R)
        speeds = unique(initCoh.speeds);
        cohs = unique(initCoh.coh);
        if strcmp('20210326a',dcp{filei}.datapath(end-8:end))
            disp(['Problems with eye trace for data associated with 20210326a; skipping...'])
        else
            passCutoff = [passCutoff; initCoh.passCutoff];
            Cnew = nan(length(binT),length(binT),length(initCoh.unitIndex),length(dir));
            for di = 1:length(dir)
                Cnew(:,:,:,di) = findCovariance(initCoh,binT,speeds,cohs,dir(di));
            end
            
            
            e = sqrt(vertcat(initCoh.eye(:).hvel).^2 + vertcat(initCoh.eye(:).vvel).^2)';
            eInit = nanmean(e(initCoh.eye_t >= initWin(1) & initCoh.eye_t <= initWin(2),:),1);
            
            % Select preferred direction
            Rtemp = [];
            Rs = [];
            Ctemp = [];
            for uniti = 1:length(initCoh.preferredDirectionRelative)
                ind = find(dir == initCoh.preferredDirectionRelative(uniti));
                Rtemp = cat(4,Rtemp,initCoh.R(:,:,:,uniti,ind));
                Ctemp = cat(3,Ctemp,Cnew(:,:,uniti,ind));
                
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        [~,condLogical] = trialSort(initCoh,dir(ind),speeds(si),NaN,cohs(ci));
                        rs = initCoh.r(:,condLogical,uniti);
                        es = eInit(condLogical);
                        Rs(:,si,ci,1,uniti) = nanmean(rs(:,es<nanmean(es)),2);
                        Rs(:,si,ci,2,uniti) = nanmean(rs(:,es>nanmean(es)),2);
                    end
                end
            end
            R = cat(4,R,Rtemp);
            Rsplit = cat(5,Rsplit,Rs);
            C = cat(3,C,Ctemp);
        
            CbehTemp = findCovarianceBehavior(initCoh,unique(initCoh.speeds),unique(initCoh.coh),0,[0 1400]);
            Cbeh = cat(3,Cbeh,CbehTemp);
            
            for uniti = 1:length(initCoh.chansIndex)
                if isempty(initCoh.location)
                    x = NaN;
                    y = NaN;
                    z = NaN;
                elseif length(initCoh.location.x)==24
                    siteIndex = chanMap(initCoh.chansIndex(uniti) == chanMap(:,1),2);
                    x = initCoh.location.x(siteIndex);
                    y = initCoh.location.y(siteIndex);
                    depth = -initCoh.location.depth(siteIndex);
                elseif length(initCoh.location.x) > 1
                    siteIndex = floor(initCoh.chansIndex(uniti)/4)+1;
                    tempIndex = find(~isnan(initCoh.location.x));
                    if siteIndex>length(tempIndex)
                        x = NaN;
                        y = NaN;
                        depth = NaN;
                    else
                        x = initCoh.location.x(tempIndex(siteIndex));
                        y = initCoh.location.y(tempIndex(siteIndex));
                        depth = -initCoh.location.depth(tempIndex(siteIndex));
                    end
                else
                    x = initCoh.location.x;
                    y = initCoh.location.y;
                    depth = -initCoh.location.depth;
                end
                locations = [locations; x,y,depth];
            end
            
            if calcCC
                [cc, shufflecc] = correlograms(initCoh,sampleRate,initCoh.unitIndex,win,shuffleN);
                for i = 1:length(initCoh.unitIndex)
                    for j = 1:length(initCoh.unitIndex)
                        CC(indx+1,j,:) = cc(:,i,j)-mean(shufflecc(:,i,j,:),4);                            % Shuffle corrected correlogram
                        shuffletemp = reshape(shufflecc(:,i,j,:),[numel(shufflecc(:,i,j,:)),1]);
                        CC_CI(indx+1,j,:) = [-1.96*std(shuffletemp),1.96*std(shuffletemp)];   % 95% confidence intervals of correlogram
                        cellID(indx+1,j,1) = filei;
                        cellID(indx+1,j,2) = initCoh.unitIndex(i);
                        cellID(indx+1,j,3) = initCoh.unitIndex(j);
                    end
                    indx = indx+1;
                end
            elseif calcRSC
                [rsc,pval] = spikeCountCorrelationWin(initCoh,'win',[150,450]);
                for i = 1:length(initCoh.unitIndex)
                    for j = 1:length(initCoh.unitIndex)
                        RSC(indx+1,j) = rsc(i,j);
                        RSC_pval(indx+1,j) = pval(i,j);
                        cellID(indx+1,j,1) = filei;
                        cellID(indx+1,j,2) = initCoh.unitIndex(i);
                        cellID(indx+1,j,3) = initCoh.unitIndex(j);
                    end
                    indx = indx+1;
                end
            else
                for i = 1:length(initCoh.unitIndex)
                    cellID(indx+1,1,1) = filei;
                    cellID(indx+1,1,2) = initCoh.unitIndex(i);
                    indx = indx+1;
                end
            end
        end
    end
end
CC = CC(1:indx,:,:);
CC_CI = CC_CI(1:indx,:,:);
cellID = cellID(1:indx,:,:);
RSC = RSC(1:indx,:);
RSC_pval = RSC_pval(1:indx,:);

%% Remove data that doesn't pass cutoff
R = R(:,:,:,logical(passCutoff));
Rsplit = Rsplit(:,:,:,:,logical(passCutoff));
C = C(:,:,logical(passCutoff));
locations = locations(logical(passCutoff),:);
cellID = cellID(logical(passCutoff),:,:);
CC = CC(logical(passCutoff),:,:);
CC_CI = CC_CI(logical(passCutoff),:,:);
RSC = RSC(logical(passCutoff),:);
RSC_pval = RSC_pval(logical(passCutoff),:);

%% Remove outlier rates
m = squeeze(max(R,[],[1,2,3,5]))*1000;
R = R(:,:,:,m<=150);
Rsplit = Rsplit(:,:,:,:,m<=150);
C = C(:,:,m<=150);
locations = locations(m<=150,:);
cellID = cellID(m<=150,:,:);
CC = CC(m<=150,:,:);
CC_CI = CC_CI(m<=150,:,:);
RSC = RSC(m<=150,:);
RSC_pval = RSC_pval(m<=150,:);

%% Mean center
R2 = R(:,:,:,:,1);
mR2 = mean(R2,[1,2,3]);
R2 = R2 - mR2;
Rsplit2 = Rsplit - repmat(permute(mR2,[1,2,3,5,4]),[1,1,1,2,1]);

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

sz = size(R);
RlowD = nan([sz(1:3),reconstructionN]);
RlowDsplit = nan([sz(1:3),2,reconstructionN]);
for si = 1:size(R,2)
    for ci = 1:size(R,3)
        RlowD(:,si,ci,:) = permute(vecs(:,1:reconstructionN)'*squeeze(R2(:,si,ci,:))',[2,3,4,1]);
        RlowDsplit(:,si,ci,1,:) = permute(vecs(:,1:reconstructionN)'*squeeze(Rsplit2(:,si,ci,1,:))',[2,3,4,5,1]);
        RlowDsplit(:,si,ci,2,:) = permute(vecs(:,1:reconstructionN)'*squeeze(Rsplit2(:,si,ci,2,:))',[2,3,4,5,1]);
    end
end

%% Use nonlinear dimensionality reduction to attempt to identify neuron 'types'
NumDimensions = 2;
Perplexity = 10;
Y = tsne(R3','Distance','cosine','NumDimensions',NumDimensions,'Perplexity',Perplexity);

switch ClusterMethod
    case 'K-means'
        NumClusters = 3;
        idx = kmeans(Y,NumClusters);

    case 'densityClust'
        dist = pdist2(Y,Y);
        percNeigh = 0.02;
        % 'Gauss' denotes the use of Gauss Kernel to compute rho, and
        % 'Cut-off' denotes the use of Cut-off Kernel.
        % For large-scale data sets, 'Cut-off' is preferable owing to computational efficiency,
        % otherwise, 'Gauss' is preferable in the case of small samples (especially with noises).
        kernel = 'Gauss';
        % set critical system parameters for DensityClust
        [dc, rho] = paraSet(dist, percNeigh, kernel);
        [NumClusters, idx, centInd, haloInd] = densityClust(dist, dc, rho, true);
end

%%
colors = colormap('lines');

%% PCA plotting
% figure
% for ri = 1:reconstructionN
%     subplot(1,reconstructionN,ri)
% %     for ci = 1:size(RlowD,3)
% %         plot(initCoh.neuron_t,RlowD(:,2,ci,ri),'Color',colors(ci,:))
% %         hold on
% %     end
%     for si = 1:size(RlowD,2)
%         plot(initCoh.neuron_t,RlowD(:,si,3,ri),'Color',colors(si,:))
%         hold on
%     end
%     xlabel('Time from motion onset (ms)')
%     ylabel(['PC ' num2str(ri)])
% end

figure
for ri = 1:reconstructionN
    for si = 1:size(RlowD,2)
        subplot(size(RlowD,2),reconstructionN,ri+(si-1)*reconstructionN)
        for ci = 1:size(RlowD,3)
            plot(initCoh.neuron_t,RlowD(:,si,ci,ri),'Color',colors(ci,:))
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

%% PCA high vs low eye speed
figure('Name','Sort by coherence')
for ri = 1:reconstructionN
    subplot(2,reconstructionN,ri)
    si = 3;
    for ci = 1:size(RlowD,3)
        plot(initCoh.neuron_t,RlowD(:,si,ci,ri),'Color',colors(ci,:))
        hold on
    end
%     xlim([0 250])
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
    subplot(2,reconstructionN,ri+reconstructionN)
    for ci = 1:size(RlowD,3)
        plot(initCoh.neuron_t,...
            RlowDsplit(:,si,ci,2,ri)-RlowDsplit(:,si,ci,1,ri),'Color',colors(ci,:))
        hold on
    end
%     ci = 3;
%     for si = 1:size(RlowD,2)
%         plot(initCoh.neuron_t,RlowDsplit(:,si,ci,2,ri)-RlowDsplit(:,si,ci,1,ri),'Color',colors(si,:))
%         hold on
% %         plot(initCoh.neuron_t,RlowDsplit(:,si,3,2,ri),'--','Color',colors(si,:))
%     end
%     xlim([0 250])
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
end

figure('Name','Sort by speed')
for ri = 1:reconstructionN
    subplot(2,reconstructionN,ri)
    ci = 3;
    for si = 1:size(RlowD,2)
        plot(initCoh.neuron_t,RlowD(:,si,ci,ri),'Color',colors(si,:))
        hold on
    end
%     xlim([0 250])
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
    subplot(2,reconstructionN,ri+reconstructionN)
    for si = 1:size(RlowD,2)
        plot(initCoh.neuron_t,RlowDsplit(:,si,ci,2,ri)-RlowDsplit(:,si,ci,1,ri),'Color',colors(si,:))
        hold on
    end
%     xlim([0 250])
    xlabel('Time from motion onset (ms)')
    ylabel(['PC ' num2str(ri)])
end

figure
ri = 1;
ind = 0;
for si = 1:size(RlowD,2)
    for ci = 1:size(RlowD,3)
        ind = ind+1;
        subplot(size(RlowD,2),size(RlowD,3),ind)
        plot(initCoh.neuron_t,sum(RlowD(:,si,ci,ri),4),'Color',[0 0 0])
        hold on
        plot(initCoh.neuron_t,sum(RlowDsplit(:,si,ci,1,ri),5),'Color',[1 0 0])
        plot(initCoh.neuron_t,sum(RlowDsplit(:,si,ci,2,ri),5),'Color',[0 0 1])
        xlim([-100 250])
        xlabel('Time from motion onset (ms)')
        ylabel(['PC ' num2str(ri)])
    end
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
subplot(2,2,1)
mC = mean(C,3);
mCbeh = nanmean(Cbeh,3);
imagesc(binT,binT,mC-diag(diag(mean(C,3))))
xlabel('Time from motion onset (ms)')
ylabel('Time from motion onset (ms)')
axis square

subplot(2,2,2)
imagesc(binT,binT,mCbeh-diag(diag(mean(Cbeh,3))))
xlabel('Time from motion onset (ms)')
ylabel('Time from motion onset (ms)')
axis square

subplot(2,2,[3 4])
binTbeh = binT(1):binT(end);
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
xlabel('Time lag (ms)')
ylabel('Relative autocorrelation')
legend({'FEFsem','Behavior'})

%% Visualize dimensionality reduction
vecInds = [1,2,3];
vColor = sqrt(vecs(:,vecInds).^2./repmat(max(vecs(:,vecInds).^2,[],2),[1,3]));

figure;
if NumDimensions > 2
    dimsTemp = [1,2,3];
    subplot(2,3,[1,2,4,5])
    scatter3(Y(:,dimsTemp(1)),Y(:,dimsTemp(2)),Y(:,dimsTemp(3)),50,vColor,'filled')
    hold on
%     scatter3(Y(:,dimsTemp(1)),Y(:,dimsTemp(2)),Y(:,dimsTemp(3)),50,colors(idx,:))
    axis square
    xlabel('tSNE 1')
    ylabel('tSNE 2')
    
    subplot(2,3,[3,6])
    for i = 1:NumClusters
        plot(initCoh.neuron_t,...
            nanmean(squeeze(R(:,3,3,idx == i))./...
            repmat(max(squeeze(R(:,3,3,idx==i)),[],1),[size(R,1) 1]),2))
        hold on
    end
    
else
    subplot(3,2,[1,2,3,4])
    scatter(Y(:,1),Y(:,2),50,colors(idx,:),'filled')
    hold on
%     scatter(Y(:,1),Y(:,2),25,vColor,'filled')
    axis square
    xlabel('tSNE 1')
    ylabel('tSNE 2')
    
    subplot(3,2,[5,6])
    for i = 1:NumClusters
        plot(initCoh.neuron_t,...
            nanmean(squeeze(R(:,3,3,idx == i))./...
            repmat(max(squeeze(R(:,3,3,idx==i)),[],1),[size(R,1) 1]),2),...
            'Color',colors(i,:))
        hold on
    end
end

figure
ind = 0;
for i = 1:NumClusters
    for si = 1:size(RlowD,2)
        ind = ind+1;
        subplot(NumClusters,size(RlowD,2),ind)
        for ci = 1:size(RlowD,3)
            tempR = nanmean(squeeze(R(:,si,ci,idx == i))./...
                repmat(max(squeeze(R(:,si,ci,idx==i)),[],1),[size(R,1) 1]),2);
            plot(initCoh.neuron_t,...
                tempR-tempR(1),...
                'Color',colors(ci,:))
            hold on
        end
        axis tight
        lims(ind,:) = axis;
    end
end
ind = 0;
for i = 1:NumClusters
    for si = 1:size(RlowD,2)
        ind = ind+1;
        subplot(NumClusters,size(RlowD,2),ind)
        axis([min(lims(:,1)),max(lims(:,2)),min(lims(:,3)),max(lims(:,4))])
    end
end

%% topography
figure;
randScale = 0.08;
vecInds = [1,2,3];
vColor = sqrt(vecs(:,vecInds).^2./repmat(max(vecs(:,vecInds).^2,[],2),[1,3]));
% randvecs = randn(size(vecs(:,vecInds)));
% vColor = sqrt(randvecs.^2./repmat(max(randvecs.^2,[],2),[1,3]));
locations2 = locations;
locations2(locations(:,1)>1,:) = [locations(locations(:,1)>1,2), locations(locations(:,1)>1,1), locations(locations(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
locationsRand = locations2 + ...
    [randScale*nanstd(locations2(:,1))*randn(size(locations,1),1), ...
    randScale*nanstd(locations2(:,2))*randn(size(locations,1),1), ...
    0*nanstd(locations2(:,3))*randn(size(locations,1),1)];                  % Add randomness to a-p and m-l locations to make 
subplot(1,2,1)
for uniti = 1:size(vColor,1)
    plot3(locationsRand(uniti,1),locationsRand(uniti,2),locationsRand(uniti,3)/1000,...
        'o','Color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',vColor(uniti,:));
    hold on
end
grid on
axis equal
xlabel('Anterior/postieror (mm)')
ylabel('Medial/lateral (mm)')
zlabel('Depth (mm)')

subplot(1,2,2)
for uniti = 1:size(vColor,1)
    plot3(locationsRand(uniti,1),locationsRand(uniti,2),locationsRand(uniti,3)/1000,...
        'o','Color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',colors(idx(uniti),:,:));
    hold on
end
grid on
axis equal
xlabel('Anterior/postieror (mm)')
ylabel('Medial/lateral (mm)')
zlabel('Depth (mm)')

%% Examine neuron 'type' conditioned cross correlations
figure
smoothWin = 1;
sigThres = 1.96;
types = [1,2];
type1 = find(idx == types(1));
CClist = [];
for uniti = 1:length(type1)
    filei = cellID(type1(uniti),1,1);
    unitIDs = cellID(type1(uniti),:,3);
    for unitj = 1:length(unitIDs)
        type2 = idx(find(cellID(:,1,1) == filei & cellID(:,1,2) ==  unitIDs(unitj)));
        if type2 == types(2)
            cc = squeeze(CC(type1(uniti),unitj,:));
            cc(win >= -1 & win <= 1) = NaN;
            cc(cc >= CC_CI(type1(uniti),unitj,1) & cc <= CC_CI(type1(uniti),unitj,2)) = NaN;
%             cc = cc - mean(cc([1:win-4 win+4:2*win]));
            CClist = [CClist cc];
        end
    end
end
mCC = nanmean(CClist,2);
% mCC(win-1:win+2) = NaN;
subplot(4,1,1)
plot(win,smooth(mCC,smoothWin))
hold on

subplot(4,1,2)
imagesc(1:size(CClist,2),win,CClist)
hold on
plotHorizontal(0);
ylim([-80 80])

types = [1,3];
type1 = find(idx == types(1));
CClist = [];
for uniti = 1:length(type1)
    filei = cellID(type1(uniti),1,1);
    unitIDs = cellID(type1(uniti),:,3);
    for unitj = 1:length(unitIDs)
        type2 = idx(find(cellID(:,1,1) == filei & cellID(:,1,2) ==  unitIDs(unitj)));
        if type2 == types(2)
            cc = squeeze(CC(type1(uniti),unitj,:));
            cc(win >= -1 & win <= 1) = NaN;
            cc(cc >= CC_CI(type1(uniti),unitj,1) & cc <= CC_CI(type1(uniti),unitj,2)) = NaN;
%             cc = cc - mean(cc([1:win-4 win+4:2*win]));
            CClist = [CClist cc];
        end
    end
end
mCC = nanmean(CClist,2);
% mCC(win-1:win+2) = NaN;

subplot(4,1,1)
plot(win,smooth(mCC,smoothWin))

subplot(4,1,3)
imagesc(1:size(CClist,2),win,CClist)
hold on
plotHorizontal(0);
ylim([-80 80])

types = [3,2];
type1 = find(idx == types(1));
CClist = [];
for uniti = 1:length(type1)
    filei = cellID(type1(uniti),1,1);
    unitIDs = cellID(type1(uniti),:,3);
    for unitj = 1:length(unitIDs)
        type2 = idx(find(cellID(:,1,1) == filei & cellID(:,1,2) ==  unitIDs(unitj)));
        if type2 == types(2)
            cc = squeeze(CC(type1(uniti),unitj,:));
            cc(win >= -1 & win <= 1) = NaN;
            cc(cc >= CC_CI(type1(uniti),unitj,1) & cc <= CC_CI(type1(uniti),unitj,2)) = NaN;
%             cc = cc - mean(cc([1:win-4 win+4:2*win]));
            CClist = [CClist cc];
        end
    end
end
mCC = nanmean(CClist,2);
% mCC(win-1:win+2) = NaN;
subplot(4,1,1)
plot(win,smooth(mCC,smoothWin))
plotHorizontal(0);
plotVertical(0);
xlim([-80 80])
% ylim([-1 1])


subplot(4,1,4)
imagesc(1:size(CClist,2),win,CClist)
hold on
plotHorizontal(0);
ylim([-80 80])

%% 

if calcCC
    % Measure connectivity
    CC2 = CC;
    for uniti = 1:size(CC_CI,1)
        for unitj = 1:size(CC_CI,2)
            CC2(uniti,unitj,CC2(uniti,unitj,:) >= CC_CI(uniti,unitj,1) & CC2(uniti,unitj,:) <= CC_CI(uniti,unitj,2)) = NaN;
        end
    end
    A = nansum(CC2(:,:,win > -50 & win < -1),3);
    B = nansum(CC2(:,:,win < 50 & win > 1),3);
    tot = nansum(CC2,3);
    deltaCC = (-A+B)./tot;
    % maxConnection = max(abs(C(:)));
    thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
    [xtemp,ytemp] = generalEllipse(1,1,'theta',thetas);
    
    figure
    circleProject = false;
    typesAll = [1 1; 1 2; 1 3; 2 1; 2 2; 2 3; 3 1; 3 2; 3 3];
    for compi = 1:size(typesAll,1)
        types = typesAll(compi,:);%[3,3];
        type1 = find(idx == types(1));
        
        tempIDs = nan(length(type1),5);
        ind = 0;
        for uniti = 1:length(type1)
            filei = cellID(type1(uniti),1,1);
            unitIDs = cellID(type1(uniti),:,3);
            for unitj = 1:length(unitIDs)
                unit2 = find(cellID(:,1,1) == filei & cellID(:,1,2) ==  unitIDs(unitj));
                type2 = idx(unit2);
                if type2 == types(2)
                    ind = ind+1;
                    if circleProject
                        tempIDs(ind,1) = xtemp(type1(uniti));%
                        tempIDs(ind,2) = xtemp(unit2);
                        tempIDs(ind,3) = ytemp(type1(uniti));%
                        tempIDs(ind,4) = ytemp(unit2);%
                        tempIDs(ind,5) = deltaCC(type1(uniti),unitj);
                    else
                        tempIDs(ind,1) = Y(type1(uniti),1);%xtemp(type1(uniti));%
                        tempIDs(ind,2) = Y(unit2,1);%xtemp(unit2);
                        tempIDs(ind,3) = Y(type1(uniti),2);%ytemp(type1(uniti));%
                        tempIDs(ind,4) = Y(unit2,2);%ytemp(unit2);%
                        tempIDs(ind,5) = deltaCC(type1(uniti),unitj);
                    end
                end
            end
        end
        tempIDs = tempIDs(1:ind,:);
        
        subplot(ceil(sqrt(size(typesAll,1))),ceil(sqrt(size(typesAll,1))),compi)
        if circleProject
            scatter(xtemp,ytemp,50,colors(idx,:),'filled')
            % scatter(xtemp,ytemp,25,vColor,'filled')
        else
            scatter(Y(:,1),Y(:,2),50,colors(idx,:),'filled')
        end
        hold on
        maxConnection = max(abs(tempIDs(:,5)));
        for pairi = 1:ind
            tempWidth = 10*sqrt(abs(tempIDs(pairi,5)))/sqrt(maxConnection);
            if tempIDs(pairi,5) > 0
                plot(tempIDs(pairi,1:2),tempIDs(pairi,3:4),...
                    'o-','Color',[0 0 0],...
                    'LineWidth',tempWidth)
            elseif tempIDs(pairi,5) < 0
                plot(tempIDs(pairi,1:2),tempIDs(pairi,3:4),...
                    'o-','Color',[1 0 0],...
                    'LineWidth',tempWidth)
            end
        end
        axis equal
        axis square
        xlabel('tSNE 1')
        ylabel('tSNE 2')
        title(num2str(typesAll(compi,:)))
    end
    
    
    vonMises = @(rho,mu,kappa)( exp(kappa*cos(rho-mu))/(2*pi*besseli(0,kappa)) );
    vonMisesKappa = (4/pi)^2;
    
    wCC = nan(size(thetas,1),NumClusters);
    wxsc = nan(NumClusters,NumClusters);
    for typei = 1:NumClusters
        thetaCen(typei) = thetas(centInd==typei);
    end
    for uniti = 1:size(thetas,1)
        filei = cellID(uniti,1,1);
        unitIDs = cellID(uniti,:,3);
        for typei = 1:NumClusters
            ccTemp = deltaCC(uniti,:);
            for unitj = 1:size(cellID,2)
                thetaInd = find(cellID(:,1,1) == filei & cellID(:,1,2) == unitIDs(unitj));
                if isempty(thetaInd)
                    thetaTemp(unitj) = NaN;
                else
                    thetaTemp(unitj) = thetas(thetaInd);
                end
            end
            wCC(uniti,typei) = nansum( vonMises(thetaTemp,thetaCen(typei),vonMisesKappa).*ccTemp );
%             wRSC(uniti,typei) = nansum(cos( (thetaTemp - thetaCen(typei))/2 ).*rscTemp)/nansum(cos( (thetaTemp - thetaCen(typei))/2 ));
        end
    end
%     sum(abs(wRSC(:))>1)
%     wRSC(abs(wRSC)>1) = NaN;
    figure
    for typei = 1:NumClusters
        subplot(1,NumClusters,typei)
        keepVec = ~isnan(wCC(:,typei)) & wCC(:,typei)~=0;
        scatter(Y(keepVec,1),Y(keepVec,2),50,wCC(keepVec,typei),'filled')
        hold on
        plot(Y(centInd==typei,1),Y(centInd==typei,2),'ko','MarkerSize',10)
        cax(typei,:) = caxis;
        axis square
    end
    for typei = 1:NumClusters
        subplot(1,NumClusters,typei)
        caxis([min(cax(:,1)) max(cax(:,2))])
        colorbar
    end
end

%%

if calcRSC
    % Measure connectivity
    RSC2 = RSC;
    RSC2(RSC_pval>0.05) = NaN;
    RSC2(RSC2 == 1) = NaN;
    thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
    [xtemp,ytemp] = generalEllipse(1,1,'theta',thetas);
    
    figure
    circleProject = false;
    [TEMP1 TEMP2] = meshgrid(1:NumClusters);
    typesAll = [TEMP1(:) TEMP2(:)];
    for compi = 1:size(typesAll,1)
        types = typesAll(compi,:);%[3,3];
        type1 = find(idx == types(1));
        
        tempIDs = nan(length(type1),5);
        ind = 0;
        for uniti = 1:length(type1)
            filei = cellID(type1(uniti),1,1);
            unitIDs = cellID(type1(uniti),:,3);
            for unitj = 1:length(unitIDs)
                unit2 = find(cellID(:,1,1) == filei & cellID(:,1,2) ==  unitIDs(unitj));
                type2 = idx(unit2);
                if type2 == types(2)
                    ind = ind+1;
                    if circleProject
                        tempIDs(ind,1) = xtemp(type1(uniti));%
                        tempIDs(ind,2) = xtemp(unit2);
                        tempIDs(ind,3) = ytemp(type1(uniti));%
                        tempIDs(ind,4) = ytemp(unit2);%
                        tempIDs(ind,5) = RSC2(type1(uniti),unitj);
                    else
                        tempIDs(ind,1) = Y(type1(uniti),1);%xtemp(type1(uniti));%
                        tempIDs(ind,2) = Y(unit2,1);%xtemp(unit2);
                        tempIDs(ind,3) = Y(type1(uniti),2);%ytemp(type1(uniti));%
                        tempIDs(ind,4) = Y(unit2,2);%ytemp(unit2);%
                        tempIDs(ind,5) = RSC2(type1(uniti),unitj);
                    end
                end
            end
        end
        tempIDs = tempIDs(1:ind,:);
        RSC_typeConditioned{compi} = tempIDs(:,5);
        mRSC(compi) = nanmean(tempIDs(:,5));
        
        subplot(ceil(sqrt(size(typesAll,1))),ceil(sqrt(size(typesAll,1))),compi)
        if circleProject
            scatter(xtemp,ytemp,50,colors(idx,:),'filled')
            % scatter(xtemp,ytemp,25,vColor,'filled')
        else
            scatter(Y(:,1),Y(:,2),50,colors(idx,:),'filled')
        end
        hold on
        maxConnection = max(abs(tempIDs(:,5)));
        for pairi = 1:ind
            tempWidth = 5*abs(tempIDs(pairi,5))/sqrt(maxConnection);
            if tempIDs(pairi,5) > 0
                plot(tempIDs(pairi,1:2),tempIDs(pairi,3:4),...
                    'o-','Color',[0 0 0],...
                    'LineWidth',tempWidth)
            elseif tempIDs(pairi,5) < 0
                plot(tempIDs(pairi,1:2),tempIDs(pairi,3:4),...
                    'o-','Color',[1 0 0],...
                    'LineWidth',tempWidth)
            end
        end
        axis equal
        axis square
        xlabel('tSNE 1')
        ylabel('tSNE 2')
        title(num2str(typesAll(compi,:)))
    end
end

%% Calculate average RSC between neurons, weighted by the cosine of the distance to center of each data type
if calcRSC
    
    vonMises = @(rho,mu,kappa)( exp(kappa*cos(rho-mu))/(2*pi*besseli(0,kappa)) );
    vonMisesKappa = (4/pi)^2;
    
    wRSC = nan(size(thetas,1),NumClusters);
    wxsc = nan(NumClusters,NumClusters);
    for typei = 1:NumClusters
        thetaCen(typei) = thetas(centInd==typei);
    end
    for uniti = 1:size(thetas,1)
        filei = cellID(uniti,1,1);
        unitIDs = cellID(uniti,:,3);
        for typei = 1:NumClusters
            rscTemp = RSC2(uniti,:);
            for unitj = 1:size(cellID,2)
                thetaInd = find(cellID(:,1,1) == filei & cellID(:,1,2) == unitIDs(unitj));
                if isempty(thetaInd)
                    thetaTemp(unitj) = NaN;
                else
                    thetaTemp(unitj) = thetas(thetaInd);
                end
            end
            wRSC(uniti,typei) = nansum( vonMises(thetaTemp,thetaCen(typei),vonMisesKappa).*rscTemp );
%             wRSC(uniti,typei) = nansum(cos( (thetaTemp - thetaCen(typei))/2 ).*rscTemp)/nansum(cos( (thetaTemp - thetaCen(typei))/2 ));
        end
    end
%     sum(abs(wRSC(:))>1)
%     wRSC(abs(wRSC)>1) = NaN;
    figure
    for typei = 1:NumClusters
        subplot(1,NumClusters,typei)
        keepVec = ~isnan(wRSC(:,typei)) & wRSC(:,typei)~=0;
        scatter(Y(keepVec,1),Y(keepVec,2),50,wRSC(keepVec,typei),'filled')
        hold on
        plot(Y(centInd==typei,1),Y(centInd==typei,2),'ko','MarkerSize',10)
        cax(typei,:) = caxis;
        axis square
    end
    for typei = 1:NumClusters
        subplot(1,NumClusters,typei)
        caxis([min(cax(:,1)) max(cax(:,2))])
        colorbar
    end
    
    figure
    subplotInd = 0;
    for typei = 1:NumClusters
        for typej = 1:NumClusters
            wxsc(typei,typej) = nansum( vonMises(thetas,thetas(centInd==typei),vonMisesKappa).*wRSC(:,typej) );
%             wxsc(typei,typej) = nansum( cos( (thetas-thetas(centInd==typei))/2 ).*wRSC(:,typej) )/nansum( cos( (thetas-thetas(centInd==typei))/2 ));
            if abs(wxsc(typei,typej))>1
                wxsc(typei,typej)=NaN;
            end
            if typej > typei
%                 subplotInd = subplotInd+1;
%                 subplot(sum(typesAll(:,1)>typesAll(:,2)),1,subplotInd)
                wwRSC = vonMises(thetas,thetas(centInd==typei),vonMisesKappa).*wRSC(:,typej);
                ntemp = histc( wwRSC.^2 ,linspace(0,0.1,50));
%                 histogram(cos((thetas -thetas(centInd==typei))/2).*wRSC(:,typej),linspace(-1,1,50))
                plot(linspace(0,0.1,50),cumsum(ntemp)/sum(ntemp))
                hold on
%                 plotVertical(wxsc(typei,typej));
                title([num2str(typei) num2str(typej)])
            end
        end
    end
    
    figure
    imagesc(wxsc)
end