function dcpBehavioralAnalysisSummary(varargin)
%% dcpBehavioralAnalysisSummary
%
%
%
%%

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'fullCalc',false)

parse(Parser,varargin{:})

fullCalc = Parser.Results.fullCalc;

if fullCalc
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
end

%% ar behavioral and physiological gain estimates 
saveFig = false;
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
saveFig = false;
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

%% ar ANOVAN 
speeds = [5, 10, 20];
cohs = [20, 60, 100];

clear init
load('/home/seth/Projects/DynamicCoherencePhysiology/ar/initCohPertResults/initCohPert20240212.mat', 'init')


% Eye speed, initiation
y = [];
spds = [];
cs = [];
SPDS = [5, 10, 20];
CS = [20, 60, 100];
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        y = [y; init.eye.init{SPDS == speeds(si), CS == cohs(ci)}(:)];
        spds = [spds; speeds(si)*ones(size(init.eye.init{SPDS == speeds(si), CS == cohs(ci)}(:)))];
        cs = [cs; cohs(ci)*ones(size(init.eye.init{SPDS == speeds(si), CS == cohs(ci)}(:)))];
    end
end
group = {spds,cs};
[summarystats.Ar.initiation.speeds.p, summarystats.Ar.initiation.speeds.tbl, summarystats.Ar.initiation.speeds.stats] = anovan(y,group,'model','interaction','varnames',{'Speed','Coherence'});


% Eye speed, sustained
y = [];
spds = [];
cs = [];
acceptVec = ismember(init.conditions.speeds,speeds) & ismember(init.conditions.coh,cohs);
y = init.eye.speed(init.t == 750,acceptVec);
spds = init.conditions.speeds(acceptVec);
cs = init.conditions.coh(acceptVec);
group = {spds,cs};
[summarystats.Ar.sustained.speeds.p, summarystats.Ar.sustained.speeds.tbl, summarystats.Ar.sustained.speeds.stats] = anovan(y,group,'model','interaction','varnames',{'Speed','Coherence'});

% Pursuit gain, initiation
pi = 2;
[summarystats.Ar.initiation.gains.Fvals, summarystats.Ar.initiation.gains.pvals, summarystats.Ar.initiation.gains.df] = approximateANOVA(init,pi);

% Pursuit gain, sustained
pi = 3;
[summarystats.Ar.sustained.gains.Fvals, summarystats.Ar.sustained.gains.pvals, summarystats.Ar.dfSustained] = approximateANOVA(init,pi);

% model comparison: linear - gain = B(1)*speed + B(2)*coh + B(3)
%                nonlinear - gain = B(1)*speed*coh + B(2)*coh + B(3)
maskConds = [0;%[10;      % speeds to test model on
             0];% 60];    % cohs to test model on
summarystats.Ar.initiation.gains.gainModel = gainModelComparison(init,2,maskConds,cohs,speeds);
summarystats.Ar.sustained.gains.gainModel = gainModelComparison(init,3,maskConds,cohs,speeds);

% Test generalization of model across time
[Cs,Ss] = meshgrid(cohs,speeds);

    % Initiation model performance predicting sustained data
    pred = summarystats.Ar.initiation.gains.gainModel.nonlinear.B(1)*Cs.*Ss + ...
        summarystats.Ar.initiation.gains.gainModel.nonlinear.B(2)*Cs + ...
        summarystats.Ar.initiation.gains.gainModel.nonlinear.B(3);
    tempGain = (init.eye.pert.res(:,:,3)-init.eye.pert.resControl(:,:,3))./(0.4.*Ss);
    summarystats.Ar.sustained.gains.gainModel.nonlinear.sseInitiationModel = ...
        sum( (pred - tempGain).^2 ,'all');
    
    % Sustained model performance predicting initiation data
    pred = summarystats.Ar.sustained.gains.gainModel.nonlinear.B(1)*Cs.*Ss + ...
        summarystats.Ar.sustained.gains.gainModel.nonlinear.B(2)*Cs + ...
        summarystats.Ar.sustained.gains.gainModel.nonlinear.B(3);
    tempGain = (init.eye.pert.res(:,:,2)-init.eye.pert.resControl(:,:,2))./(0.4.*Ss);
    summarystats.Ar.initiation.gains.gainModel.nonlinear.sseSustainedModel = ...
        sum( (pred - tempGain).^2 ,'all');

% Specific test of difference between high speed, high coherence gain
% during initiatoin is smaller than high speed, high coherence gain during
% sustained pursuit
signal(1) = (init.eye.pert.res(3,3,2)-init.eye.pert.resControl(3,3,2))./(0.4.*20);
signal(2) = (init.eye.pert.res(3,3,3)-init.eye.pert.resControl(3,3,3))./(0.4.*20);
STE(1) = sqrt(init.eye.pert.resSTE(3,3,2).^2 + init.eye.pert.resControlSTE(3,3,2).^2)/(0.4*20);
STE(2) = sqrt(init.eye.pert.resSTE(3,3,3).^2 + init.eye.pert.resControlSTE(3,3,3).^2)/(0.4*20);
summarystats.Ar.SpecificTests.highspeed.t = (signal(2) - signal(1))/sqrt(STE(2)^2 + STE(1)^2);
DoF = length(init.eye.init{3,3})/2 -1;
summarystats.Ar.SpecificTests.highspeed.pval = 1-tcdf(summarystats.Ar.SpecificTests.highspeed.t,DoF);    % Original hypothesis is one tailed (sustained > than intiation)
summarystats.Ar.SpecificTests.highspeed.DoF = DoF;

% Specific test of difference between low speed, high coherence gain
% during initiatoin is smaller than low speed, high coherence gain during
% sustained pursuit
signal(1) = (init.eye.pert.res(1,3,2)-init.eye.pert.resControl(1,3,2))./(0.4.*5);
signal(2) = (init.eye.pert.res(1,3,3)-init.eye.pert.resControl(1,3,3))./(0.4.*5);
STE(1) = sqrt(init.eye.pert.resSTE(1,3,2).^2 + init.eye.pert.resControlSTE(1,3,2).^2)/(0.4*5);
STE(2) = sqrt(init.eye.pert.resSTE(1,3,3).^2 + init.eye.pert.resControlSTE(1,3,3).^2)/(0.4*5);
summarystats.Ar.SpecificTests.lowspeed.t = (signal(2) - signal(1))/sqrt(STE(2)^2 + STE(1)^2);
DoF = length(init.eye.init{1,3})/2 -1;
summarystats.Ar.SpecificTests.lowspeed.pval(1) = 1-tcdf(summarystats.Ar.SpecificTests.lowspeed.t,DoF);    % Original hypothesis is one tailed (sustained > than intiation)
summarystats.Ar.SpecificTests.lowspeed.pval(2) = 2*(1-tcdf(abs(summarystats.Ar.SpecificTests.lowspeed.t),DoF));    % Post-hoc test if there is any difference
summarystats.Ar.SpecificTests.lowspeed.DoF = DoF;  
    
%% fr ANOVAN 
speeds = [5, 10, 20];
cohs = [20, 60, 100];
clear init
load('/home/seth/Projects/DynamicCoherencePhysiology/fr/initCohPertResults/initCohPert20240212.mat', 'init')


% Eye speed, initiation
y = [];
spds = [];
cs = [];
SPDS = [5, 10, 20];
CS = [20, 60, 100];
for si = 1:length(speeds)
    for ci = 1:length(cohs)
        y = [y; init.eye.init{SPDS == speeds(si), CS == cohs(ci)}(:)];
        spds = [spds; speeds(si)*ones(size(init.eye.init{SPDS == speeds(si), CS == cohs(ci)}(:)))];
        cs = [cs; cohs(ci)*ones(size(init.eye.init{SPDS == speeds(si), CS == cohs(ci)}(:)))];
    end
end
group = {spds,cs};
[summarystats.Fr.initiation.speeds.p, summarystats.Fr.initiation.speeds.tbl, summarystats.Fr.initiation.speeds.stats] = anovan(y,group,'model','interaction','varnames',{'Speed','Coherence'});


% Eye speed, sustained
y = [];
spds = [];
cs = [];
acceptVec = ismember(init.conditions.speeds,speeds) & ismember(init.conditions.coh,cohs);
y = init.eye.speed(init.t == 750,acceptVec);
spds = init.conditions.speeds(acceptVec);
cs = init.conditions.coh(acceptVec);
group = {spds,cs};
[summarystats.Fr.sustained.speeds.p, summarystats.Fr.sustained.speeds.tbl, summarystats.Fr.sustained.speeds.stats] = anovan(y,group,'model','interaction','varnames',{'Speed','Coherence'});

% Pursuit gain, initiation
pi = 2;
[summarystats.Fr.initiation.gains.Fvals, summarystats.Fr.initiation.gains.pvals, summarystats.Fr.initiation.gains.df] = approximateANOVA(init,pi);

% Pursuit gain, sustained
pi = 3;
[summarystats.Fr.sustained.gains.Fvals, summarystats.Fr.sustained.gains.pvals, summarystats.Fr.dfSustained] = approximateANOVA(init,pi);

% model comparison: linear - gain = B(1)*speed + B(2)*coh + B(3)
%                nonlinear - gain = B(1)*speed*coh + B(2)*coh + B(3)
maskConds = [0;%[10;      % speeds to test model on
             0];% 60];    % cohs to test model on
summarystats.Fr.initiation.gains.gainModel = gainModelComparison(init,2,maskConds,cohs,speeds);
summarystats.Fr.sustained.gains.gainModel = gainModelComparison(init,3,maskConds,cohs,speeds);

% Test generalization of model across time
[Cs,Ss] = meshgrid(cohs,speeds);

    % Initiation model performance predicting sustained data
    pred = summarystats.Fr.initiation.gains.gainModel.nonlinear.B(1)*Cs.*Ss + ...
        summarystats.Fr.initiation.gains.gainModel.nonlinear.B(2)*Cs + ...
        summarystats.Fr.initiation.gains.gainModel.nonlinear.B(3);
    tempGain = (init.eye.pert.res(:,:,3)-init.eye.pert.resControl(:,:,3))./(0.4.*Ss);
    summarystats.Fr.sustained.gains.gainModel.nonlinear.sseInitiationModel = ...
        sum( (pred - tempGain).^2 ,'all');
    
    % Sustained model performance predicting initiation data
    pred = summarystats.Fr.sustained.gains.gainModel.nonlinear.B(1)*Cs.*Ss + ...
        summarystats.Fr.sustained.gains.gainModel.nonlinear.B(2)*Cs + ...
        summarystats.Fr.sustained.gains.gainModel.nonlinear.B(3);
    tempGain = (init.eye.pert.res(:,:,2)-init.eye.pert.resControl(:,:,2))./(0.4.*Ss);
    summarystats.Fr.initiation.gains.gainModel.nonlinear.sseSustainedModel = ...
        sum( (pred - tempGain).^2 ,'all');

% Specific test of difference between high speed, high coherence gain
% during initiatoin is smaller than high speed, high coherence gain during
% sustained pursuit
signal(1) = (init.eye.pert.res(3,3,2)-init.eye.pert.resControl(3,3,2))./(0.4.*20);
signal(2) = (init.eye.pert.res(3,3,3)-init.eye.pert.resControl(3,3,3))./(0.4.*20);
STE(1) = sqrt(init.eye.pert.resSTE(3,3,2).^2 + init.eye.pert.resControlSTE(3,3,2).^2)/(0.4*20);
STE(2) = sqrt(init.eye.pert.resSTE(3,3,3).^2 + init.eye.pert.resControlSTE(3,3,3).^2)/(0.4*20);
summarystats.Fr.SpecificTests.highspeed.t = (signal(2) - signal(1))/sqrt(STE(2)^2 + STE(1)^2);
DoF = length(init.eye.init{3,3})/2 -1;
summarystats.Fr.SpecificTests.highspeed.pval = 1-tcdf(summarystats.Fr.SpecificTests.highspeed.t,DoF);    % Original hypothesis is one tailed (sustained > than intiation)
summarystats.Fr.SpecificTests.highspeed.DoF = DoF;

% Specific test of difference between low speed, high coherence gain
% during initiatoin is smaller than low speed, high coherence gain during
% sustained pursuit
signal(1) = (init.eye.pert.res(1,3,2)-init.eye.pert.resControl(1,3,2))./(0.4.*5);
signal(2) = (init.eye.pert.res(1,3,3)-init.eye.pert.resControl(1,3,3))./(0.4.*5);
STE(1) = sqrt(init.eye.pert.resSTE(1,3,2).^2 + init.eye.pert.resControlSTE(1,3,2).^2)/(0.4*5);
STE(2) = sqrt(init.eye.pert.resSTE(1,3,3).^2 + init.eye.pert.resControlSTE(1,3,3).^2)/(0.4*5);
summarystats.Fr.SpecificTests.lowspeed.t = (signal(2) - signal(1))/sqrt(STE(2)^2 + STE(1)^2);
DoF = length(init.eye.init{1,3})/2 -1;
summarystats.Fr.SpecificTests.lowspeed.pval(1) = 1-tcdf(summarystats.Fr.SpecificTests.lowspeed.t,DoF);    % Original hypothesis is one tailed (sustained > than intiation)
summarystats.Fr.SpecificTests.lowspeed.pval(2) = 2*(1-tcdf(abs(summarystats.Fr.SpecificTests.lowspeed.t),DoF));    % Post-hoc test if there is any difference
summarystats.Fr.SpecificTests.lowspeed.DoF = DoF;
    
%% Plot summary of model performance across animals
figure('Name','Behavioral model comparison')
subplot(1,2,1)
plot(sqrt(summarystats.Ar.initiation.gains.gainModel.linear.sse),...
    sqrt(summarystats.Ar.initiation.gains.gainModel.nonlinear.sse),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Ar linear vs nonlinear gain model, initiation')
hold on
plot(sqrt(summarystats.Fr.initiation.gains.gainModel.linear.sse),...
    sqrt(summarystats.Fr.initiation.gains.gainModel.nonlinear.sse),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Fr linear vs nonlinear gain model, initiation')
plot(sqrt(summarystats.Ar.sustained.gains.gainModel.linear.sse),...
    sqrt(summarystats.Ar.sustained.gains.gainModel.nonlinear.sse),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Ar linear vs nonlinear gain model, sustained')
plot(sqrt(summarystats.Fr.sustained.gains.gainModel.linear.sse),...
    sqrt(summarystats.Fr.sustained.gains.gainModel.nonlinear.sse),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Fr linear vs nonlinear gain model, sustained')
xlabel('Linear SSE')
ylabel('Nonlinear SSE')
title('Linear vs nonlinear model performance')
plotUnity;
axis square;

subplot(1,2,2)
plot(sqrt(summarystats.Ar.initiation.gains.gainModel.nonlinear.sse),...
    sqrt(summarystats.Ar.initiation.gains.gainModel.nonlinear.sseSustainedModel),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Ar initiation vs sustained model, initiation')
hold on
plot(sqrt(summarystats.Fr.initiation.gains.gainModel.nonlinear.sse),...
    sqrt(summarystats.Fr.initiation.gains.gainModel.nonlinear.sseSustainedModel),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Fr initiation vs sustained model, initiation')
plot(sqrt(summarystats.Ar.sustained.gains.gainModel.nonlinear.sse),...
    sqrt(summarystats.Ar.sustained.gains.gainModel.nonlinear.sseInitiationModel),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Ar sustained vs initiation model, sustained')
plot(sqrt(summarystats.Fr.sustained.gains.gainModel.nonlinear.sse),...
    sqrt(summarystats.Fr.sustained.gains.gainModel.nonlinear.sseInitiationModel),...
    'o','Color',[0 0 0],'MarkerFaceColor',[0 0 0],...
    'DisplayName','Fr sustained vs initiation model, sustained')
xlabel('Dynamic SSE')
ylabel('Static SSE')
title('Linear vs nonlinear model performance')
plotUnity;
axis square

%% ar initialCoh

%% fr initialCoh

%% ar dynamicCoh
detectPoorPursuit.On = true;
detectPoorPursuit.threshold = 1.5;
[dyn, gain] = dynamicCohBehavioralAnalysis('ar','dcpObjectsFile','dcpObjects20210406.mat');

%% fr dynamicCoh


%% Functions

%% approximateANOVA
function [Fvals, pvals, df] = approximateANOVA(init,pi)
    means = init.eye.pert.res(:,:,pi)-init.eye.pert.resControl(:,:,pi);
    std_devs = sqrt(init.eye.pert.resSTE(:,:,pi).^2.*init.eye.pert.N(:,:,pi) +...
        init.eye.pert.resControlSTE(:,:,pi).^2.*init.eye.pert.N(:,:,pi));
    sample_sizes = init.eye.pert.N(:,:,pi);

    % Compute Grand Mean
    grand_mean = mean(means(:));

    % Compute Sum of Squares (SS) for Main Effects and Interaction
    SS_A = sum(sum(sample_sizes .* (mean(means, 2) - grand_mean).^2));
    SS_B = sum(sum(sample_sizes .* (mean(means, 1) - grand_mean).^2));
    SS_AB = sum(sum(sample_sizes .* ((means - grand_mean).^2))) - (SS_A + SS_B);

    % Compute Degrees of Freedom
    df.Speed = size(means, 1) - 1;  % df for Factor A (3 levels → 2 df)
    df.Coherence = size(means, 2) - 1;  % df for Factor B (3 levels → 2 df)
    df.SpeedXCoherence = df.Speed * df.Coherence;        % df for Interaction (2 × 2 = 4)
    df.error = sum(sample_sizes(:)) - numel(means); % Approximate error df

    % Compute Mean Squares
    MS_A = SS_A / df.Speed;
    MS_B = SS_B / df.Coherence;
    MS_AB = SS_AB / df.SpeedXCoherence;

    % Approximate F-tests (assuming homogeneity of variance)
    pooled_variance = mean(std_devs(:).^2); % Estimate error variance
    Fvals.Speed = MS_A / pooled_variance;
    Fvals.Coherence = MS_B / pooled_variance;
    Fvals.SpeedXCoherence = MS_AB / pooled_variance;
    
    pvals.Speed = 1-fcdf(Fvals.Speed,df.Speed,df.error);
    pvals.Coherence = 1-fcdf(Fvals.Coherence,df.Coherence,df.error);
    pvals.SpeedXCoherence = 1-fcdf(Fvals.SpeedXCoherence,df.SpeedXCoherence,df.error);
    
%% Compare regression models of gain
function gainModel = gainModelComparison(init,pi,maskConds,cohs,speeds)
    [Cs,Ss] = meshgrid(cohs,speeds);
    Mask = ~(ismember(Cs,maskConds(2,:)) & ismember(Ss,maskConds(1,:)));
    tempRes = (init.eye.pert.res(:,:,pi)-init.eye.pert.resControl(:,:,pi))./(0.4.*Ss);
    tempSTE = sqrt(init.eye.pert.resSTE(:,:,pi).^2 + init.eye.pert.resControlSTE(:,:,pi).^2)./(0.4*Ss);

    [gainModel.linear.B, gainModel.linear.BINT, ~,~, gainModel.linear.STATS] = regress(tempRes(Mask(:)),[Ss(Mask(:)) Cs(Mask(:)) ones(size(Cs(Mask(:))))]);
    [gainModel.nonlinear.B, gainModel.nonlinear.BINT, ~,~, gainModel.nonlinear.STATS] = regress(tempRes(Mask(:)),[Ss(Mask(:)).*Cs(Mask(:)) Cs(Mask(:)) ones(size(Cs(Mask(:))))]);

    if any(~Mask(:))
        gainModel.linear.sse = sum( ([Ss(~Mask(:)) Cs(~Mask(:)) ones(size(Cs(~Mask(:))))]*gainModel.linear.B - tempRes(~Mask(:))).^2 );
        gainModel.nonlinear.sse = sum( ([Ss(~Mask(:)).*Cs(~Mask(:)) Cs(~Mask(:)) ones(size(Cs(~Mask(:))))]*gainModel.nonlinear.B - tempRes(~Mask(:))).^2 );
    else
        gainModel.linear.sse = sum( ([Ss(:) Cs(:) ones(size(Cs(:)))]*gainModel.linear.B - tempRes(:)).^2 );
        gainModel.nonlinear.sse = sum( ([Ss(:).*Cs(:) Cs(:) ones(size(Cs(:)))]*gainModel.nonlinear.B - tempRes(:)).^2 );
    end
