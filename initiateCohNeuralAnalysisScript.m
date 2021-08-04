%% Defaults
binT = 0:100:1400;
dir = [0 180];

%% Get mean and covariance of each unit
load('dcpObjects20210406.mat');
sourceDirectory = '/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle';

C = [];
passCutoff = [];
R = [];
for filei = 1:length(dcp)
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    if ~isempty(initCoh.R)
        
        passCutoff = [passCutoff; initCoh.passCutoff];
        Cnew = findCovariance(initCoh,binT,speeds,cohs,0);
        
        C = cat(3,C,Cnew);
        
        % Select preferred direction
        Rtemp = [];
        for uniti = 1:length(initCoh.preferredDirectionRelative)
            ind = find(dir == initCoh.preferredDirectionRelative(uniti));
            Rtemp = cat(4,Rtemp,initCoh.R(:,:,:,:,ind));
        end
        R = cat(4,R,Rtemp);
        
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


%%
colors = colormap('lines');
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
mC = mean(C,3)-diag(diag(mean(C,3)));
imagesc(binT,binT,mC)
axis square

subplot(1,2,2)
for bi = 1:length(binT)
    plot(binT(bi+1:end),mC(bi,bi+1:end))
    hold on
end