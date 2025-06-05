function [weightedAverage, standardAverage, interpolatedWeightedAverage, interpolatedStandardAverage, conditions, t, figureHandles] = analyzeMTcoherencePertubationData(varargin)
%% analyzeMTcoherencePerturbationData
%
%   Analyzing MT data from Behling and Lisberger (2023). Data is converted
%   from .jld2 to .mat by ~/Projects/DyanmicCoherencePhysiology/ar/mtdata/jld2ToMAT.jl
%
%       To get default data, navigate to ~/Projects/DyanmicCoherencePhysiology/ar/mtdata/
%       in the terminal and type
%           julia jld2ToMAT MTdataProcessed.jld2 /mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/MTdata/MTdataProcessed.mat
%
%
%%

%% Defaults
conditions_default = [8 100   0 0 1;          % Sp = 8, Coh = 100, cohPulse = 0, speedPulse = 0, prefD = True
                      8  30   0 0 1;          % Sp = 8, Coh = 30, cohPulse = 0, speedPulse = 0, prefD = True
                      8  10   0 0 1;          % Sp = 8, Coh = 10, cohPulse = 0, speedPulse = 0, prefD = True
                      8  10  90 0 1;          % Sp = 8, Coh = 10, cohPulse = 90, speedPulse = 0, prefD = True
                      8 100 -90 0 1;          % Sp = 8, Coh = 10, cohPulse = -90, speedPulse = 0, prefD = True
                      8  10  20 0 1;          % Sp = 8, Coh = 10, cohPulse = 20, speedPulse = 0, prefD = True
                      8 100 -70 0 1];         % Sp = 8, Coh = 100, cohPulse = -70, speedPulse = 0, prefD = True

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'MTdataFile','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/MTdata/MTdataProcessed.mat')
addParameter(Parser,'t',-400:5:900)
addParameter(Parser,'initialCohs',[100, 30, 10, 10, 100, 10, 100])
addParameter(Parser,'pulseAmplitudes',[0, 0, 0, 90, -90, 20, -70])
addParameter(Parser,'speed',[8, 16])
addParameter(Parser,'speedPulseData',false)
addParameter(Parser,'preferredDirection',1)
addParameter(Parser,'transitionTime',540)
addParameter(Parser,'timeconstant',20)
addParameter(Parser,'removeBaseline',true)
addParameter(Parser,'smoothData',true)
addParameter(Parser,'snames',{'Ar','Di'})
addParameter(Parser,'plotFlg',true)

parse(Parser,varargin{:})

MTdataFile = Parser.Results.MTdataFile;
t = Parser.Results.t;
initialCohs = Parser.Results.initialCohs;
pulseAmplitudes = Parser.Results.pulseAmplitudes;
speed = Parser.Results.speed;
speedPulseData = Parser.Results.speedPulseData;
preferredDirection = Parser.Results.preferredDirection;
transitionTime = Parser.Results.transitionTime;
timeconstant = Parser.Results.timeconstant;
removeBaseline = Parser.Results.removeBaseline;
smoothData = Parser.Results.smoothData;
snames = Parser.Results.snames;
plotFlg = Parser.Results.plotFlg;

transitionInd = find(t==transitionTime);


%% Load MT data
load(MTdataFile)

%% Create conditions table

% Conditons table should be of the form:
% [Speed, initialCoh, pulseAmplitude, speedPulse, prefD]

controlCohs = unique([initialCohs + pulseAmplitudes]);

ind = 0;
for si = 1:length(speed)
    for condi = 1:length(initialCohs)
        ind = ind+1;
        conditions(ind,:) = [speed(si), initialCohs(condi), pulseAmplitudes(condi), speedPulseData, preferredDirection];
    end
end

% Lia = ismember(controlCohs,conditions(pulseAmplitudes==0,2));
% controlsNotIncluded = controlCohs(~Lia);
% for condi = 1:length(controlsNotIncluded)
%     conditions(condi+ind,:) = [speed(si), controlsNotIncluded(condi), 0, speedPulseData, preferredDirection];
% end


%% Convert to standard arrays
spdall = vertcat(spdall{:});
cohPulse = vertcat(cohPulse{:});
speedPulse = vertcat(speedPulse{:});
psall = vertcat(psall{:});


%% Analyze MT responses to coherence changes
for mti = 1:length(mdall)
    cellName = mdall{mti};
    pulseData(mti,1) = any(strcmp(mdall,cellName) & cohPulse > 0);
    dataFR(:,mti) = dataFRall{mti}(1:length(t));
    nameAccept(mti,1) = any(strcmp(mdall{mti}(1:2),snames));
end

for condi = 1:size(conditions,1)
    logicalInds(condi,:) = spdall == conditions(condi,1) & ...
        cohall == conditions(condi,2) & ...
        cohPulse == conditions(condi,3) & ...
        speedPulse == conditions(condi,4) & ...
        prefD == conditions(condi,5) & ...
        pulseData & ...
        nameAccept;
end

% Now average normalized firing rates across neurons for the relavant conditions
% sp8_coh100_noPulse = spdall == 8 & ...
%     cohall == 100 & ...
%     cohPulse == 0 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;
% 
% sp8_coh30_noPulse = spdall == 8 & ...
%     cohall == 30 & ...
%     cohPulse == 0 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;
% 
% sp8_coh10_noPulse = spdall == 8 & ...
%     cohall == 10 & ...
%     cohPulse == 0 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;
% 
% sp8_coh10_Pulse = spdall == 8 & ...
%     cohall == 10 & ...
%     cohPulse == 90 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;
% 
% sp8_coh100_Pulse = spdall == 8 & ...
%     cohall == 100 & ...
%     cohPulse == -90 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;
% 
% sp8_coh10_Pulse2 = spdall == 8 & ...
%     cohall == 10 & ...
%     cohPulse == 20 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;
% 
% sp8_coh100_Pulse2 = spdall == 8 & ...
%     cohall == 100 & ...
%     cohPulse == -70 & ...
%     speedPulse == 0 & ...
%     prefD == 1 & ...
%     pulseData;

%% Compute averege and weighted averages
for condi = 1:size(conditions,1)
    weightedAverage(:,condi) = sum(dataFR(:,logicalInds(condi,:))./max(dataFR(:,logicalInds(condi,:))) .* ...
        (log2(psall(logicalInds(condi,:))) - log2(10))',2)/sum(logicalInds(condi,:));
    
    standardAverage(:,condi) = mean(dataFR(:,logicalInds(condi,:))./max(dataFR(:,logicalInds(condi,:))),2);
    
    if removeBaseline
        weightedAverage(:,condi) = weightedAverage(:,condi) - mean(weightedAverage(t<0,condi),'all');
        standardAverage(:,condi) = standardAverage(:,condi) - mean(standardAverage(t<0,condi),'all');
    end
   
end

%% Expected response based on interpolation from original coherence to perturbed coherence
interpolatedWeightedAverage = nan(size(weightedAverage));

condsToInterpolate = find(conditions(:,3));

for condi = 1:length(condsToInterpolate)
    origStream = weightedAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
        conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0);
    finalStream = weightedAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
        conditions(:,2) == conditions(condsToInterpolate(condi),2)+conditions(condsToInterpolate(condi),3) & ...
        conditions(:,3) == 0);
    interpolatedWeightedAverage(:,condsToInterpolate(condi)) = frInterpolate(origStream,...
        finalStream,transitionInd,0);
    
    origStream = standardAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
        conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0);
    finalStream = standardAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
        conditions(:,2) == conditions(condsToInterpolate(condi),2)+conditions(condsToInterpolate(condi),3) & ...
        conditions(:,3) == 0);
    interpolatedStandardAverage(:,condsToInterpolate(condi)) = frInterpolate(origStream,...
        finalStream,transitionInd,0);    
end


 
if smoothData
    for condi = 1:size(conditions,1)
        weightedAverage(:,condi) = lowpass(weightedAverage(:,condi),1/timeconstant,1000);
        standardAverage(:,condi) = lowpass(standardAverage(:,condi),1/timeconstant,1000);
        if ~any(isnan(interpolatedWeightedAverage(:,condi)))
            interpolatedWeightedAverage(:,condi) = lowpass(interpolatedWeightedAverage(:,condi),1/timeconstant,1000);
        end
        if ~any(isnan(interpolatedStandardAverage(:,condi)))
            interpolatedStandardAverage(:,condi) = lowpass(interpolatedStandardAverage(:,condi),1/timeconstant,1000);
        end
    end
end

%% Normalize by variance before motion onset
wa_pow = std(weightedAverage(t<0,:),[],'all');
sa_pow = std(standardAverage(t<0,:),[],'all');
weightedAverage = weightedAverage / wa_pow;
interpolatedWeightedAverage = interpolatedWeightedAverage / wa_pow;
standardAverage = standardAverage / sa_pow;
interpolatedStandardAverage = interpolatedStandardAverage / sa_pow;

%% Normalize by maximum across conditions
wa_max = max([weightedAverage; standardAverage],[],'all');
sa_max = max([weightedAverage; standardAverage],[],'all');
weightedAverage = weightedAverage ./ wa_max;
interpolatedWeightedAverage = interpolatedWeightedAverage ./ wa_max;
standardAverage = standardAverage ./ sa_max;
interpolatedStandardAverage = interpolatedStandardAverage ./ sa_max;

%% Ploting
if plotFlg
    
    %% Plot each coherence perturbation condition data, corresponding data without perturbations and interpolated data
    figureHandles{1} = figure;
    for condi = 1:length(condsToInterpolate)
        
        subplot(length(condsToInterpolate),2,1+(condi-1)*2)
        origStream = weightedAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0);
        finalStream = weightedAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2)+conditions(condsToInterpolate(condi),3) & ...
            conditions(:,3) == 0);
        plot(t,origStream,...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0,2))])
        hold on
        plot(t,finalStream,...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2)+conditions(condsToInterpolate(condi),3) & ...
            conditions(:,3) == 0,2))])
        plot(t,interpolatedWeightedAverage(:,condsToInterpolate(condi)),...
            'DisplayName','Interpolation')
        plot(t,weightedAverage(:,condsToInterpolate(condi)),...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0,2)) ...
            ', cohPulse = ' num2str(conditions(condsToInterpolate(condi),3))])
        plotVertical(500);
        
        
        subplot(length(condsToInterpolate),2,2+(condi-1)*2)
        origStream = standardAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0);
        finalStream = standardAverage(:,conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2)+conditions(condsToInterpolate(condi),3) & ...
            conditions(:,3) == 0);
        plot(t,origStream,...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0,2))])
        hold on
        plot(t,finalStream,...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2)+conditions(condsToInterpolate(condi),3) & ...
            conditions(:,3) == 0,2))])
        plot(t,interpolatedStandardAverage(:,condsToInterpolate(condi)),...
            'DisplayName','Interpolation')
        plot(t,standardAverage(:,condsToInterpolate(condi)),...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0,2)) ...
            ', cohPulse = ' num2str(conditions(condsToInterpolate(condi),3))])
        plotVertical(500);
    end
    
    %%
    figureHandles{2} = figure;
    for condi = 1:length(condsToInterpolate)
        
        subplot(1,2,1)
        plot(t,weightedAverage(:,condsToInterpolate(condi))-interpolatedWeightedAverage(:,condsToInterpolate(condi)),...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0,2)) ...
            ', cohPulse = ' num2str(conditions(condsToInterpolate(condi),3))])
        hold on
        
        subplot(1,2,2)
        plot(t,standardAverage(:,condsToInterpolate(condi))-interpolatedStandardAverage(:,condsToInterpolate(condi)),...
            'DisplayName',['Inital coherence = ' num2str(conditions(conditions(:,1) == conditions(condsToInterpolate(condi),1) & ...
            conditions(:,2) == conditions(condsToInterpolate(condi),2) & conditions(:,3) == 0,2)) ...
            ', cohPulse = ' num2str(conditions(condsToInterpolate(condi),3))])
        hold on
    end
    subplot(1,2,1)
    plotVertical(500);
    plotHorizontal(0);
    xlabel('Time from motion onset (ms)')
    ylabel('Effect not captured by interpolation')
    
    subplot(1,2,2)
    plotVertical(500);
    plotHorizontal(0);
    xlabel('Time from motion onset (ms)')
    ylabel('Effect not captured by interpolation')
    
end

%% Functions
function interpolatedFR = frInterpolate(stream1,stream2,transitionInd,timeconstant)
    %%
    interpolatedFR = stream1;
    streamDiff = stream2(transitionInd:end)-stream1(transitionInd:end);
    interpolatedFR(transitionInd:end) = interpolatedFR(transitionInd:end) + ...
        streamDiff;
    if timeconstant > 0
        interpolatedFR = lowpass(interpolatedFR,1/timeconstant,1000);
    end