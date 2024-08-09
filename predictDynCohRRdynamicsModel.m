function predictDynCohRRdynamicsModel(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'neuronTypingFile',NaN)
addParameter(Parser,'fitReducedRankDynamicsFile',NaN)
addParameter(Parser,'MTtoFEFregressionFile',NaN)
addParameter(Parser,'speedInd',2)
addParameter(Parser,'cohInds',[2 1 3])
addParameter(Parser,'transitionTime',490)
addParameter(Parser,'timeconstant',20)
addParameter(Parser,'inputGain',[8 2])

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
neuronTypingFile = Parser.Results.neuronTypingFile;
fitReducedRankDynamicsFile = Parser.Results.fitReducedRankDynamicsFile;
MTtoFEFregressionFile = Parser.Results.MTtoFEFregressionFile;
speedInd = Parser.Results.speedInd;
cohInds = Parser.Results.cohInds;
transitionTime = Parser.Results.transitionTime;
timeconstant = Parser.Results.timeconstant;
inputGain = Parser.Results.inputGain;


%% Check for specific files as input, otherwise set default files
if any(isnan(neuronTypingFile))
    switch subject
        case 'ar'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
        case 'fr'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';            
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(fitReducedRankDynamicsFile))
    switch subject
        case 'ar'
%             fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240613.mat';
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240729.mat';
        case 'fr'
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/fitReducedRankModel20240613.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(MTtoFEFregressionFile))
    switch subject
        case 'ar'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/ar/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        case 'fr'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/fr/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

%% Load the data
neuronTyping = load(neuronTypingFile,'subjects','Y','cellID','gainRegression','idx','Rinit','NumClusters','initCoh');
rrModel = load(fitReducedRankDynamicsFile,'modelFEF','theoreticalInput','kappas');
MTtoFEF = load(MTtoFEFregressionFile,'vOpt');
t = rrModel.modelFEF.t;
modelFEF = rrModel.modelFEF;

%% Create weighted MT response for 60% coh to transitions to 20 and 100% coh
inputInitCoh = rrModel.theoreticalInput;
input = inputInitCoh(:,:,speedInd,cohInds(1));
transitionInd = find(t == transitionTime);
for i = 2:length(cohInds)
    input(1,:,i) = inputStreamInterpolate(inputInitCoh(1,:,speedInd,cohInds(1)),...
        inputInitCoh(1,:,speedInd,cohInds(i)),...
        transitionInd,timeconstant,[1,1]);
    input(2,:,i) = inputStreamInterpolate(inputInitCoh(2,:,speedInd,cohInds(1)),...
        inputInitCoh(2,:,speedInd,cohInds(i)),...
        transitionInd,timeconstant,inputGain);
end
% input(1,:,:) = inputInitCoh(1,:,speedInd,cohInds);
input(3,:,:) = inputInitCoh(3,:,3,1:3);

%% Run simulation of DynCoh
for i = 1:length(cohInds)
    [~,kappas(:,:,i)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
        't',modelFEF.t,'us',input(:,:,i),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
end

%% Plot
c = lines;
figure
for ki = 1:size(kappas,1)
    subplot(1,size(kappas,1),ki)
    for ci = 1:size(rrModel.kappas,4)
        plot(t,rrModel.kappas(ki,:,speedInd,ci),'k')
        hold on
    end
    for ci = 2:length(cohInds)
        plot(t,kappas(ki,:,ci),'Color',c(ci,:))
    end
    xlabel('Time from motion onset (ms)')
end

figure
for inputi = 1:size(input,1)
    subplot(1,size(input,1),inputi)
    plot(t,squeeze(inputInitCoh(inputi,:,speedInd,cohInds)),'k')
    hold on
    ind = 1:size(inputInitCoh,3);
%     plot(t,squeeze(inputInitCoh(inputi,:,ind~=speedInd,3)),'k--')
%     plot(t,squeeze(inputInitCoh(inputi,:,ind~=speedInd,2)),'k--')
%     plot(t,squeeze(inputInitCoh(inputi,:,ind~=speedInd,1)),'k--')
    for ci = 2:length(cohInds)
        plot(t,squeeze(input(inputi,:,ci)),'Color',c(ci,:))
    end
end

%% Analyze MT responses to coherence changes


for mti = 1:length(mdall)
    cellName = mdall{mti};
    pulseData(mti,1) = any(strcmp(mdall,cellName) & cohPulse > 0);
end


% Now average normalized firing rates across neurons for the relavant conditions
sp8_coh100_noPulse = spdall == 8 & ...
    cohall == 100 & ...
    cohPulse == 0 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

sp8_coh30_noPulse = spdall == 8 & ...
    cohall == 30 & ...
    cohPulse == 0 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

sp8_coh10_noPulse = spdall == 8 & ...
    cohall == 10 & ...
    cohPulse == 0 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

sp8_coh10_Pulse = spdall == 8 & ...
    cohall == 10 & ...
    cohPulse == 90 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

sp8_coh100_Pulse = spdall == 8 & ...
    cohall == 100 & ...
    cohPulse == -90 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

sp8_coh10_Pulse2 = spdall == 8 & ...
    cohall == 10 & ...
    cohPulse == 20 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

sp8_coh100_Pulse2 = spdall == 8 & ...
    cohall == 100 & ...
    cohPulse == -70 & ...
    speedPulse == 0 & ...
    prefD == 1 & ...
    pulseData;

%%

figure
% Weighted by log2 pref speed
subplot(2,2,1)
weightedSum = sum(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)).*(log2(psall(sp8_coh100_noPulse)) - log2(10))',2)/sum(sp8_coh100_noPulse);
% plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
%     'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 100, No pulse'])
hold on
weightedSum = sum(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)).*(log2(psall(sp8_coh10_noPulse)) - log2(10))',2)/sum(sp8_coh10_noPulse);
plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
    'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 10, No pulse'])
weightedSum = sum(dataFRall(:,sp8_coh10_Pulse)./max(dataFRall(:,sp8_coh10_Pulse)).*(log2(psall(sp8_coh10_Pulse)) - log2(10))',2)/sum(sp8_coh10_Pulse);
plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
    'Color',c(3,:),'DisplayName',['Sp=8, Coh = 10, +90 pulse'])
weightedSum = sum(dataFRall(:,sp8_coh100_Pulse)./max(dataFRall(:,sp8_coh100_Pulse)).*(log2(psall(sp8_coh100_Pulse)) - log2(10))',2)/sum(sp8_coh100_Pulse);
% plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
%     'Color',c(2,:),'DisplayName',['Sp=8, Coh = 100, -90 pulse'])
axis tight
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')
ylim([-0.5 0.5])
plotVertical(500);
title('Effect of 90% Coh pulse')

subplot(2,2,3)
% weightedSum = sum(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)).*(log2(psall(sp8_coh100_noPulse)) - log2(10))',2)/sum(sp8_coh100_noPulse);
% plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
%     'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 100, No pulse'])
hold on
weightedSum = sum(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)).*(log2(psall(sp8_coh10_noPulse)) - log2(10))',2)/sum(sp8_coh10_noPulse);
plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
    'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 10, No pulse'])
weightedSum = sum(dataFRall(:,sp8_coh10_Pulse2)./max(dataFRall(:,sp8_coh10_Pulse2)).*(log2(psall(sp8_coh10_Pulse2)) - log2(10))',2)/sum(sp8_coh10_Pulse2);
plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
    'Color',c(3,:),'DisplayName',['Sp=8, Coh = 10, +20 pulse'])
% weightedSum = sum(dataFRall(:,sp8_coh100_Pulse2)./max(dataFRall(:,sp8_coh100_Pulse2)).*(log2(psall(sp8_coh100_Pulse2)) - log2(10))',2)/sum(sp8_coh100_Pulse2);
% plot(mt_t,weightedSum - mean(weightedSum(mt_t<0)),...
%     'Color',c(2,:),'DisplayName',['Sp=8, Coh = 100, -70 pulse'])
axis tight
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')
ylim([-0.5 0.5])
plotVertical(500);
title('Effect of 20%/70% Coh pulse')

% Equal weighting
subplot(2,2,2)
% plot(mt_t,mean(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),2) - ...
%     mean(dataFRall(mt_t<0,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),'all'),...
%     'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 100, No pulse'])
hold on
plot(mt_t,mean(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),'all'),...
    'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 10, No pulse'])
plot(mt_t,mean(dataFRall(:,sp8_coh10_Pulse)./max(dataFRall(:,sp8_coh10_Pulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh10_Pulse)./max(dataFRall(:,sp8_coh10_Pulse)),'all'),...
    'Color',c(3,:),'DisplayName',['Sp=8, Coh = 10, +90 pulse'])
% plot(mt_t,mean(dataFRall(:,sp8_coh100_Pulse)./max(dataFRall(:,sp8_coh100_Pulse)),2) - ...
%     mean(dataFRall(mt_t<0,sp8_coh100_Pulse)./max(dataFRall(:,sp8_coh100_Pulse)),'all'),...
%     'Color',c(2,:),'DisplayName',['Sp=8, Coh = 100, -90 pulse'])
axis tight
ylim([-0.5 0.5])
plotVertical(500);
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')
title('Effect of 90% Coh pulse')

subplot(2,2,4)
% plot(mt_t,mean(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),2) - ...
%     mean(dataFRall(mt_t<0,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),'all'),...
%     'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 100, No pulse'])
hold on
plot(mt_t,mean(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),'all'),...
    'Color',[0 0 0],'DisplayName',['Sp=8, Coh = 10, No pulse'])
plot(mt_t,mean(dataFRall(:,sp8_coh10_Pulse2)./max(dataFRall(:,sp8_coh10_Pulse2)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh10_Pulse2)./max(dataFRall(:,sp8_coh10_Pulse2)),'all'),...
    'Color',c(3,:),'DisplayName',['Sp=8, Coh = 10, +20 pulse'])
% plot(mt_t,mean(dataFRall(:,sp8_coh100_Pulse2)./max(dataFRall(:,sp8_coh100_Pulse2)),2) - ...
%     mean(dataFRall(mt_t<0,sp8_coh100_Pulse2)./max(dataFRall(:,sp8_coh100_Pulse2)),'all'),...
%     'Color',c(2,:),'DisplayName',['Sp=8, Coh = 100, -70 pulse'])
axis tight
ylim([-0.5 0.5])
plotVertical(500);
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')
title('Effect of 20%/70% Coh pulse')

%%
figure
subplot(2,2,1)
weightedSum = sum(dataFRall(:,sp8_coh10_Pulse)./max(dataFRall(:,sp8_coh10_Pulse)).*(log2(psall(sp8_coh10_Pulse)) - log2(10))',2)/sum(sp8_coh10_Pulse);
weightedSum = weightedSum - mean(weightedSum(mt_t<0));
weightedSum2 = sum(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)).*(log2(psall(sp8_coh10_noPulse)) - log2(10))',2)/sum(sp8_coh10_noPulse);
weightedSum2 = weightedSum2 - mean(weightedSum2(mt_t<0));
plot(mt_t,(weightedSum-weightedSum2)*0.5,...
    'DisplayName',['(Sp=8, Coh = 100, -90 pulse) minus (Sp=8, Coh = 100, noPulse)'])
hold on
weightedSum = sum(dataFRall(:,sp8_coh100_Pulse)./max(dataFRall(:,sp8_coh100_Pulse)).*(log2(psall(sp8_coh100_Pulse)) - log2(10))',2)/sum(sp8_coh100_Pulse);
weightedSum = weightedSum - mean(weightedSum(mt_t<0));
weightedSum2 = sum(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)).*(log2(psall(sp8_coh100_noPulse)) - log2(10))',2)/sum(sp8_coh100_noPulse);
weightedSum2 = weightedSum2 - mean(weightedSum2(mt_t<0));
plot(mt_t,(weightedSum-weightedSum2)*0.5,...
    'DisplayName',['(Sp=8, Coh = 100, -90 pulse) minus (Sp=8, Coh = 10, noPulse)'])
axis tight
ylim([-0.25 0.25])
plotHorizontal(0);
plotVertical(500);
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')

subplot(2,2,2)
hold on
plot(mt_t,mean(dataFRall(:,sp8_coh10_Pulse)./max(dataFRall(:,sp8_coh10_Pulse)),2) - ...
    mean(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh10_Pulse)./max(dataFRall(:,sp8_coh10_Pulse)),'all') + ...
    mean(dataFRall(mt_t<0,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),'all'),...
    'DisplayName',['(Sp=8, Coh = 10, +90 pulse) minus (Sp=8, Coh = 10, noPulse)'])
plot(mt_t,mean(dataFRall(:,sp8_coh100_Pulse)./max(dataFRall(:,sp8_coh100_Pulse)),2) - ...
    mean(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh100_Pulse)./max(dataFRall(:,sp8_coh100_Pulse)),'all') + ...
    mean(dataFRall(mt_t<0,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),'all'),...
    'DisplayName',['(Sp=8, Coh = 100, -90 pulse) minus (Sp=8, Coh = 100, noPulse)'])
axis tight
ylim([-0.25 0.25])
plotHorizontal(0);
plotVertical(500);
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')

subplot(2,2,3)
weightedSum = sum(dataFRall(:,sp8_coh10_Pulse2)./max(dataFRall(:,sp8_coh10_Pulse2)).*(log2(psall(sp8_coh10_Pulse2)) - log2(10))',2)/sum(sp8_coh10_Pulse2);
weightedSum = weightedSum - mean(weightedSum(mt_t<0));
weightedSum2 = sum(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)).*(log2(psall(sp8_coh10_noPulse)) - log2(10))',2)/sum(sp8_coh10_noPulse);
weightedSum2 = weightedSum2 - mean(weightedSum2(mt_t<0));
plot(mt_t,(weightedSum-weightedSum2)*0.5,...
    'DisplayName',['(Sp=8, Coh = 100, 20% pulse) minus (Sp=8, Coh = 100, noPulse)'])
hold on
weightedSum = sum(dataFRall(:,sp8_coh100_Pulse2)./max(dataFRall(:,sp8_coh100_Pulse2)).*(log2(psall(sp8_coh100_Pulse2)) - log2(10))',2)/sum(sp8_coh100_Pulse2);
weightedSum = weightedSum - mean(weightedSum(mt_t<0));
weightedSum2 = sum(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)).*(log2(psall(sp8_coh100_noPulse)) - log2(10))',2)/sum(sp8_coh100_noPulse);
weightedSum2 = weightedSum2 - mean(weightedSum2(mt_t<0));
plot(mt_t,(weightedSum-weightedSum2)*0.5,...
    'DisplayName',['(Sp=8, Coh = 100, -70% pulse) minus (Sp=8, Coh = 10, noPulse)'])
axis tight
ylim([-0.25 0.25])
plotHorizontal(0);
plotVertical(500);
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')

subplot(2,2,4)
plot(mt_t,mean(dataFRall(:,sp8_coh10_Pulse2)./max(dataFRall(:,sp8_coh10_Pulse2)),2) - ...
    mean(dataFRall(:,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh10_Pulse2)./max(dataFRall(:,sp8_coh10_Pulse2)),'all') + ...
    mean(dataFRall(mt_t<0,sp8_coh10_noPulse)./max(dataFRall(:,sp8_coh10_noPulse)),'all'),...
    'DisplayName',['(Sp=8, Coh = 10, +20 pulse) minus (Sp=8, Coh = 10, noPulse)'])
hold on
plot(mt_t,mean(dataFRall(:,sp8_coh100_Pulse2)./max(dataFRall(:,sp8_coh100_Pulse2)),2) - ...
    mean(dataFRall(:,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),2) - ...
    mean(dataFRall(mt_t<0,sp8_coh100_Pulse2)./max(dataFRall(:,sp8_coh100_Pulse2)),'all') + ...
    mean(dataFRall(mt_t<0,sp8_coh100_noPulse)./max(dataFRall(:,sp8_coh100_noPulse)),'all'),...
    'DisplayName',['(Sp=8, Coh = 100, -70 pulse) minus (Sp=8, Coh = 100, noPulse)'])
axis tight
ylim([-0.25 0.25])
plotHorizontal(0);
plotVertical(500);
xlabel('Time from motion onset (ms)')
ylabel('Summed normalized responses')




%% Functions
function input = inputStreamInterpolate(stream1,stream2,transitionInd,timeconstant,inputGain)
    %%
    t = 0:100;
    f = exp(-t/timeconstant);
    f = f/sum(f);
    input = stream1;
    streamDiff = stream2(:,transitionInd:end)-stream1(:,transitionInd:end);
    for inputi = 1:size(input,1)
        streamDiff_pos = zeros(size(streamDiff(inputi,:)));
        streamDiff_pos(streamDiff(inputi,:) > 0) = streamDiff(inputi,streamDiff(inputi,:) > 0);
        streamDiff_neg = zeros(size(streamDiff(inputi,:)));
        streamDiff_neg(streamDiff(inputi,:) < 0) = streamDiff(inputi,streamDiff(inputi,:) < 0);
        input(inputi,transitionInd:end) = input(inputi,transitionInd:end) + ...
            inputGain(1)*streamDiff_pos + inputGain(2)*streamDiff_neg;
        %     input(:,transitionInd-length(t):end) = filter(f,1,input(:,transitionInd-length(t):end),[],2);
        input(inputi,:) = lowpass(input(inputi,:),1/timeconstant,1000);
    end