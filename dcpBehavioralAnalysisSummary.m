%% dcpBehavioralAnalysisSummary
%
%
%
%%

%% ar initialCoh pertrubation
detectPoorPursuit.On = true;
detectPoorPursuit.threshold = 1.5;
[init,gain,gainModel] = initialCohPertBehavioralAnalysis('ar','dcpObjectsFile','dcpObjectsPertTemp',...
    'detectPoorPursuit',detectPoorPursuit,'saveResults',false,'saveFigures',false);

%% fr initialCoh perturbation
detectPoorPursuit.On = true;
detectPoorPursuit.threshold = 1.5;
[init,gain,gainModel] = initialCohPertBehavioralAnalysis('fr','dcpObjectsFile','dcpObjectsPert20230929',...
    'detectPoorPursuit',detectPoorPursuit,'saveResults',false,'saveFigures',false);

%% ar behavioral and physiological gain estimates 
saveFig = true;
speeds = [5, 10, 20];
cohs = [20, 60, 100];

load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat', 'init')
slips = -squeeze(init.eye.mean(init.t == 750,:,:)) + speeds';

for pi = 1:3
    for ci = 1:length(cohs)
        gain(:,ci,pi) = (init.eye.pert.res(:,ci,pi)'-init.eye.pert.resControl(:,ci,pi)')./(0.4*speeds);
        gain95CI(:,ci,pi) = sqrt(init.eye.pert.resSTE(:,ci,pi)'.^2 + init.eye.pert.resControlSTE(:,ci,pi)'.^2)./(0.4*speeds)*1.64;
    end
    
    if pi == 2
        % Initiation
        gPhys(:,:,pi) = gain(:,:,pi).*speeds'*0.4./log2(1.4);
        gPhys95CI(:,:,pi) = gain95CI(:,:,pi).*speeds'*0.4./log2(1.4);
    elseif pi == 3
        % Closed loop
        gPhys(:,:,pi) = gain(:,:,pi).*speeds'*0.4./log2(1+0.4*speeds'./slips);
        gPhys95CI(:,:,pi) = gain95CI(:,:,pi).*speeds'*0.4./log2(1+0.4*speeds'./slips);
    end
end

h = figure('Name','Estimate of physiological gain');
subplot(2,2,1)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gain(:,ci,2),gain95CI(:,ci,2),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 50 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Behavioral gain')
set(gca,'TickDir','out')

subplot(2,2,2)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gain(:,ci,3),gain95CI(:,ci,3),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 600 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Behavioral gain')
set(gca,'TickDir','out')

subplot(2,2,3)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gPhys(:,ci,2),gPhys95CI(:,ci,2),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 50 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Physiological gain')
set(gca,'TickDir','out')

subplot(2,2,4)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gPhys(:,ci,3),gPhys95CI(:,ci,3),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 600 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Physiological gain')
set(gca,'TickDir','out')

if saveFig
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/ar/initCohBehavior/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    savefig(h,[saveLocation '/physiologicalGainByTargetSpeed'])
end

%% fr behavioral and physiological gain estimates 
saveFig = true;
speeds = [5, 10, 20];
cohs = [20, 60, 100];

load('/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212.mat', 'init')
slips = -squeeze(init.eye.mean(init.t == 750,:,:)) + speeds';

for pi = 1:3
    for ci = 1:length(cohs)
        gain(:,ci,pi) = (init.eye.pert.res(:,ci,pi)'-init.eye.pert.resControl(:,ci,pi)')./(0.4*speeds);
        gain95CI(:,ci,pi) = sqrt(init.eye.pert.resSTE(:,ci,pi)'.^2 + init.eye.pert.resControlSTE(:,ci,pi)'.^2)./(0.4*speeds)*1.64;
    end
    
    if pi == 2
        % Initiation
        gPhys(:,:,pi) = gain(:,:,pi).*speeds'*0.4./log2(1.4);
        gPhys95CI(:,:,pi) = gain95CI(:,:,pi).*speeds'*0.4./log2(1.4);
    elseif pi == 3
        % Closed loop
        gPhys(:,:,pi) = gain(:,:,pi).*speeds'*0.4./log2(1+0.4*speeds'./slips);
        gPhys95CI(:,:,pi) = gain95CI(:,:,pi).*speeds'*0.4./log2(1+0.4*speeds'./slips);
    end
end

h = figure('Name','Estimate of physiological gain');
subplot(2,2,1)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gain(:,ci,2),gain95CI(:,ci,2),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 50 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Behavioral gain')
set(gca,'TickDir','out')

subplot(2,2,2)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gain(:,ci,3),gain95CI(:,ci,3),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 600 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Behavioral gain')
set(gca,'TickDir','out')

subplot(2,2,3)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gPhys(:,ci,2),gPhys95CI(:,ci,2),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 50 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Physiological gain')
set(gca,'TickDir','out')

subplot(2,2,4)
for ci = 1:size(gPhys,2)
    errorbar(speeds,gPhys(:,ci,3),gPhys95CI(:,ci,3),[],'o-',...
        'Color',1-cohs(ci)/100*ones(1,3),'MarkerFaceColor',1-cohs(ci)/100*ones(1,3),...
        'DisplayName',[num2str(cohs(ci)) '%'])
    hold on
end
title(['Perturbation time = 600 ms'])
xlim([0.8*min(speeds) 1.2*max(speeds)])
xlabel('Target speed (deg/s)')
ylabel('Physiological gain')
set(gca,'TickDir','out')

if saveFig
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/fr/initCohBehavior/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    savefig(h,[saveLocation '/physiologicalGainByTargetSpeed'])
end

%% ar initialCoh

%% fr initialCoh

%% ar dynamicCoh
detectPoorPursuit.On = true;
detectPoorPursuit.threshold = 1.5;
[dyn, gain] = dynamicCohBehavioralAnalysis('ar','dcpObjectsFile','dcpObjects20210406.mat');

%% fr dynamicCoh