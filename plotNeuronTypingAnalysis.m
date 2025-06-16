function plotNeuronTypingAnalysis(subject,varargin)
%%
%
%   Loads resultts of neuronTypingAnalysis_streamlined and plots the
%   results.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'date','20240530')
addParameter(Parser,'saveFigures',false)

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
date = Parser.Results.date;
saveFigures = Parser.Results.saveFigures;

%% Load results
fileNameLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/'...
    subject '/neuronTypingAnalysis/neuronTyping' date '.mat'];

load(fileNameLocation)

%% Plot
%% Example neurons
    exUnits = [77,140;...
        67,186;...
        69, 81;...
        70,156;...
        75, 80;...
        57,350];
    for uniti = 1:size(exUnits,1)
        
        exIndex = exUnits(uniti,:);
        cellID2 = squeeze(cellID(:,1,1:2));
        listIndex = find(ismember(cellID2, exIndex, 'rows'));
        hExUnits(uniti) = figure('Name',['Unit ' num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)}],...
            'Position',[1553 976 973 341]);
        ylims = [Inf,-Inf];
        if initCohCollate
            for speedi = 1:length(speeds)
                subplot(2,length(speeds),speedi)
                for ci = 1:length(cohs)
                    plot(initCoh.neuron_t,Rinit(:,speedi,ci,listIndex)*1000,'Color',initColors(ci,:),...
                        'DisplayName',['Neuron ' num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)} ', Speed = ' num2str(speeds(speedi)) ' deg/s, Coh = ' num2str(cohs(ci)) '%']);
                    hold on
                end
                xlabel('Time from motion onset (ms)')
                ylabel('Spikes/s')
                title('initCoh')
                axis tight
                tempLims = ylim;
                ylims(1) = min([ylims(1),tempLims(1)]);
                ylims(2) = max([ylims(2),tempLims(2)]);
            end
        end
        
        subplot(2,length(speeds),2+length(speeds))
        if dynCohCollate
            for seqi = 1:5
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),Rdyn(:,seqi,listIndex)*1000,'Color',colors(seqi,:),...
                    'DisplayName',['Neuron ' num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)} ', Sequence = ' num2str(seqi)]);
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
        title('dynCoh')
        axis tight
        tempLims = ylim;
        ylims(1) = min([ylims(1),tempLims(1)]);
        ylims(2) = max([ylims(2),tempLims(2)]);
        
        newax = [-100 1350 ylims];
        for speedi = 1:length(speeds)
            subplot(2,length(speeds),speedi)
            axis(newax)
        end
        subplot(2,length(speeds),2+length(speeds))
        axis(newax)
        text((newax(2)-newax(1))*.5+newax(1),(newax(4)-newax(3))*.1+newax(3),...
            [num2str(exIndex) ', subject ' subjects{cellID(listIndex,1,4)} ])
    end
    
    %% PSTH averaged accross neurons in a cluster
    for i = 1:NumClusters
        if initCohCollate
            a = Rinit(initCoh.neuron_t<=1350,:,:,idx == i);
        else
            a = nan(size(Rinit(:,:,:,idx == i)));
        end
        b = nan(size(Rdyn(:,:,idx == i)));
        c = nan(size(Rdyn(:,:,idx == i)));
        if dynCohCollate
            for seqi = 1:size(Rdyn,2)
                t20 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 20);
                if ~isempty(t20)
                    b(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,:) = ...
                        Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),seqi,idx == i) - ...
                        Rdyn(dynCoh.neuron_t>=min(t20) & dynCoh.neuron_t<=max(t20),5,idx == i);
                end
                
                t100 = dynCoh.eye_t(dynCoh.coh(:,seqi) == 100);
                if ~isempty(t100)
                    c(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,:) = ...
                        Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),seqi,idx == i) - ...
                        Rdyn(dynCoh.neuron_t>=min(t100) & dynCoh.neuron_t<=max(t100),5,idx == i);
                end
            end
        end
%         
%         figure('Name',['Cluster ' num2str(i)])
%         subplot(2,1,1)
%         for seqi = 1:5
%             plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(c(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
%             hold on
%             plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(b(:,seqi,:)),2)*1000,'Color',colors(seqi,:))
%         end
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)-a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(1,:))
%         hold on
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)-a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
%         plotHorizontal(0);
%         plotVertical([450 750 1050]);
%         xlabel('Time from motion onset (ms)')
%         ylabel('Excess spike/s')
%         
%         subplot(2,1,2)
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
%         hold on
%         plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(:,5,idx == i)),2)*1000,'Color',colors(5,:))
%         plotVertical([450 750 1050]);
%         xlabel('Time from motion onset (ms)')
%         ylabel('Spikes/s')
        
        hClusterPSTH(i) = figure('Name',['Cluster ' num2str(i)]);
        subplot(3,1,1)
        if dynCohCollate
            for seqi = 1:5
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(dynCoh.neuron_t<=1350,seqi,idx == i)),2)*1000,...
                    'Color',colors(seqi,:))
                hold on
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),nanmean(squeeze(Rdyn(dynCoh.neuron_t<=1350,seqi,idx == i)),2)*1000,...
                    'Color',colors(seqi,:))
            end
        end
        plotVertical([450 750 1050]);
        xlabel('Time from motion onset (ms)')
        ylabel('Spike/s')
        
        subplot(3,1,2)
        if initCohCollate
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,1,:)),2)*1000,'Color',initColors(1,:))
            hold on
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,2,:)),2)*1000,'Color',initColors(2,:))
            plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,speeds == initSpeed,3,:)),2)*1000,'Color',initColors(3,:))
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Spike/s')
        
        subplot(3,1,3)
        if initCohCollate
            for si = 1:length(speeds)
                plot(initCoh.neuron_t(initCoh.neuron_t<=1350),nanmean(squeeze(a(:,si,2,:)),2)*1000,'Color',speedColors(si,:))
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Spikes/s')
    end
       
    %% Gain vs average response in cluster
    for i = 1:NumClusters
        gvrh(i) = figure('Name',['Cluster ' num2str(i)]);
        for subjecti = 1:length(subjects)
            subplot(2,length(subjects),subjecti)
            if dynCohCollate
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn{subjecti}(eye_tDyn{subjecti} >= 450 & eye_tDyn{subjecti} <= 450,:,:),1));
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                        hold on
                    end
                end
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn{subjecti}(eye_tDyn{subjecti} >= 750 & eye_tDyn{subjecti} <= 750,:,:),1));
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                    end
                end
                eyeSpeedTemp = squeeze(nanmean(meanEyeSpeedDyn{subjecti}(eye_tDyn{subjecti} >= 1050 & eye_tDyn{subjecti} <= 1050,:,:),1));
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(eyeSpeedTemp(seqi),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                    end
                end
            end
            if initCohCollate
                for si = 1:length(speeds)
                    eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed{subjecti}(eye_t{subjecti} >= 750 & eye_t{subjecti} <= 750,:,:),1));
                    initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                    plot(eyeSpeedTemp(si,:),initRatesTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 750 ms'])
                    hold on
                end
                if includeEarlyInitCohPertTime
                    for si = 1:length(speeds)
                        eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed{subjecti}(eye_t{subjecti} >= 150 & eye_t{subjecti} <= 150,:,:),1));
                        initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                        plot(eyeSpeedTemp(si,:),initRatesTemp,...
                            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 50 ms'])
                    end
                end
                
            end
            xlabel('Eye speed (deg/s)')
            ylabel('Spikes/s')
            title(subjects(subjecti))
            set(gca,'TickDir','out')
            
            subplot(2,length(subjects),subjecti + length(subjects))
            if dynCohCollate
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 400 & dynCoh.neuron_t <= 500,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(gain{subjecti}(seqi,2),dynRatesTemp(seqi),'d','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                        hold on
                    end
                end
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 700 & dynCoh.neuron_t <= 800,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(gain{subjecti}(seqi,3),dynRatesTemp(seqi),'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                    end
                end
                dynRatesTemp = nanmean(squeeze(nanmean(Rdyn(dynCoh.neuron_t >= 1000 & dynCoh.neuron_t <= 1100,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                if length(dynRatesTemp) == length(sequences)
                    for seqi = 1:length(sequences)
                        plot(gain{subjecti}(seqi,4),dynRatesTemp(seqi),'s','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                            'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                    end
                end
            end
            if initCohCollate
                for si = 1:length(speeds)
                    initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 700 & initCoh.neuron_t <= 800,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                    plot(squeeze(initGain{subjecti}(si,:,3)),initRatesTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 750 ms'])
                    hold on
                end
                if includeEarlyInitCohPertTime
                    for si = 1:length(speeds)
                        initRatesTemp = nanmean(squeeze(nanmean(Rinit(initCoh.neuron_t >= 100 & initCoh.neuron_t <= 200,si,:,idx'==i & cellID(:,1,4)==subjecti),1)),2)*1000;
                        plot(squeeze(initGain{subjecti}(si,:,2)),initRatesTemp,...
                            'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ', 50 ms'])
                    end
                end
                
            end
            xlabel('Behavioral gain (unitless)')
            ylabel('Spikes/s')
            title(subjects(subjecti))
            set(gca,'TickDir','out')
        end
    end
    
    %% Gain decoding
    hGainDecodingFromClusters = figure('Name','Gain decoding','Position',[1633 927 888 395]);
    for subjecti = 1:length(subjects)
        subplot(length(subjects),2,1+2*(subjecti-1))
        if initCohCollate && dynCohCollate
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi),gainRegression(subjecti).yhat(seqi),...
                    'd','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+length(sequences)),gainRegression(subjecti).yhat(seqi+length(sequences)),...
                    'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+2*length(sequences)),gainRegression(subjecti).yhat(seqi+2*length(sequences)),...
                    's','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                hold on
            end
            for si = 1:length(speeds)
                for ci = 1:length(cohs)
                    plot(gainRegression(subjecti).y(3*length(sequences)+(si-1)*length(cohs)+ci),...
                        gainRegression(subjecti).yhat(3*length(sequences)+(si-1)*length(cohs)+ci),...
                        'o','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 750 ms'])
                    hold on
                end
            end
            if includeEarlyInitCohPertTime
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        plot(gainRegression(subjecti).y(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+ci),...
                            gainRegression(subjecti).yhat(3*length(sequences)+3*length(cohs)+(si-1)*length(cohs)+ci),...
                            'd','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 50 ms'])
                        hold on
                    end
                end
            end
        elseif initCohCollate
            for si = 1:length(speeds)
                for ci = 1:length(cohs)
                    plot(gainRegression(subjecti).y((si-1)*length(cohs)+ci),...
                        gainRegression(subjecti).yhat((si-1)*length(cohs)+ci),...
                        'o','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                        'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 750 ms'])
                    hold on
                end
            end
            if includeEarlyInitCohPertTime
                for si = 1:length(speeds)
                    for ci = 1:length(cohs)
                        plot(gainRegression(subjecti).y(3*length(cohs)+(si-1)*length(cohs)+ci),...
                            gainRegression(subjecti).yhat(3*length(cohs)+(si-1)*length(cohs)+ci),...
                            'd','Color',initColors(ci,:),'MarkerFaceColor',initColors(ci,:),...
                            'DisplayName',[subjects{subjecti} ', initCoh speed = ' num2str(speeds(si)) ' deg/s, coh = ' num2str(cohs(ci)), ' 50 ms'])
                        hold on
                    end
                end
            end
        elseif dynCohCollate
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi),gainRegression(subjecti).yhat(seqi),...
                    'd','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 450 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+length(sequences)),gainRegression(subjecti).yhat(seqi+length(sequences)),...
                    'o','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 750 ms'])
                hold on
            end
            for seqi = 1:length(sequences)
                plot(gainRegression(subjecti).y(seqi+2*length(sequences)),gainRegression(subjecti).yhat(seqi+2*length(sequences)),...
                    's','Color',colors(seqi,:),'MarkerFaceColor',colors(seqi,:),...
                    'DisplayName',[subjects{subjecti} ', dynCoh seq ' num2str(seqi) ', 1050 ms'])
                hold on
            end           
        end
        plotUnity;
        axis equal
        axis tight
        xlabel('Behaioral gain')
        ylabel('Estimated gain')
        set(gca,'TickDir','out')
        
        subplot(length(subjects),2,2+2*(subjecti-1))
        stem(regressorClusters,gainRegression(subjecti).B(1:end-1))
        xlabel('Cluster')
        ylabel('Regression weight')
        axis square
        set(gca,'TickDir','out')
    end
    
    %% Topology of gain
    hGainTopology = figure('Name','Gain topology','Position',[1030 669 1312 468]);
    for subjecti = 1:length(subjects)
        subplot(1,length(subjects),subjecti)
        tempGain = nan(size(Y,1),1);
        for typei = 1:NumClusters
            tempGain(idx'==typei) = gainRegression(subjecti).B(typei)/(max(gainRegression(subjecti).B(1:NumClusters)) - min(gainRegression(subjecti).B(1:NumClusters)));
        end
        scatter(Y(:,1),Y(:,2),abs(tempGain)*200,tempGain,'filled')
        hold on
        xlabel('tSNE 1')
        ylabel('tSNE 2')
        axis equal
        axis square
        set(gca,'TickDir','out')
        colormap autumn
        title(subjects{subjecti})
    end
        
        
            
    
    %% Visualize dimensionality reduction
    
    hFunctionalTopology = figure('Name','Functional topography','Position',[1396 220 560 1109]);
    subjectMarkers = {'o','square','diamond','^','<','>','pentagram','hexagram'};
    subplot(4,2,[1,2,3,4])
    for typei = 1:NumClusters
        for subjecti = 1:length(subjects)
            plot(Y(idx'==typei & cellID(:,1,4) == subjecti,1),Y(idx'==typei & cellID(:,1,4) == subjecti,2),...
                subjectMarkers{subjecti},'Color',colorWheel(typei,:),...
                'MarkerFaceColor',colorWheel(typei,:),...
                'DisplayName',[subjects{subjecti} ', cluster ' num2str(typei)]);
            hold on
        end
    end
    axis square
    xlabel('tSNE 1')
    ylabel('tSNE 2')
    
    subplot(4,2,[5,6])
    if initCohCollate
        for i = 1:NumClusters
            plot(initCoh.neuron_t,...
                nanmean(squeeze(Rinit(:,3,3,idx == i))./...
                repmat(max(squeeze(Rinit(:,3,3,idx==i)),[],1),[size(Rinit,1) 1]),2),...
                'Color',colorWheel(i,:),...
                'DisplayName',['InitCoh, speed = 20 deg/s, Coh  = 100%, cluster ' num2str(i)])
            hold on
        end
    end
    title('Response by cluster in initCoh, 20 deg/s target and 100% coherence')
    xlabel('Time from motion onset (ms)')
    ylabel('Normalized response')
    
    subplot(4,2,[7,8])
    if dynCohCollate
        for i = 1:NumClusters
            plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
                nanmean(squeeze(Rdyn(:,5,idx == i))./...
                repmat(max(squeeze(Rdyn(:,5,idx==i)),[],1),[size(Rdyn,1) 1]),2),...
                'Color',colorWheel(i,:),...
                'DisplayName',['DynCoh, speed = 10 deg/s, sequence 5, cluster ' num2str(i)])
            hold on
        end
    end
    title('Response by cluster in dynCoh, control sequence')
    xlabel('Time from motion onset (ms)')
    ylabel('Normalized response')
    
    %% Responses for each speed and coherence in initCoh
    hClusterPSTHsInitCoh = figure('Name','Coherence effect, initCoh','Position',[19 196 3*570 1133]);
    ind = 0;
    for i = 1:NumClusters
        for si = 1:size(Rinit,2)
            ind = ind+1;
            subplot(NumClusters,size(Rinit,2),ind)
            if initCohCollate
                for ci = 1:size(Rinit,3)
                    tempR = nanmean(squeeze(Rinit(:,si,ci,idx == i)),2);
                    plot(initCoh.neuron_t,...
                        tempR,...
                        'Color',initColors(ci,:),...
                        'DisplayName',['Speed = ' num2str(speeds(si)) ' deg/s, Coh = ' num2str(cohs(ci)) '%'])
                    hold on
                end
            end
            axis tight
            lims(ind,:) = axis;
        end
    end
    ind = 0;
    for i = 1:NumClusters
        for si = 1:size(Rinit,2)
            ind = ind+1;
            subplot(NumClusters,size(Rinit,2),ind)
            axis([min(lims(:,1)),max(lims(:,2)),min(lims(:,3)),max(lims(:,4))])
            xlabel('Time from motion onset (ms)')
            ylabel('Normalized response')
        end
    end
    
    %% Responses for each cluster in dynCoh
    hClusterPSTHsDynCoh = figure('Name','Coherence effect, dynCoh','Position',[1956 196 570 1133]);
    for i = 1:NumClusters
        subplot(NumClusters,1,i)
        if dynCohCollate
            for seqi = 1:size(Rdyn,2)
                tempR = nanmean(squeeze(Rdyn(:,seqi,idx == i))./...
                    repmat(max(squeeze(Rdyn(:,seqi,idx==i)),[],1),[size(Rdyn,1) 1]),2);
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),...
                    tempR,...
                    'Color',colors(seqi,:),'DisplayName',['Sequence ' num2str(seqi)])
                hold on
            end
        end
        axis tight
        lims2(i,:) = axis;
    end
    for i = 1:NumClusters
        subplot(NumClusters,1,i)
        for seqi = 1:size(Rdyn,2)
            axis([min(lims2(:,1)),max(lims2(:,2)),min(lims2(:,3)),max(lims2(:,4))])
            plotVertical([450 750 1050]);
            xlabel('Time from motion onset (ms)')
            ylabel('Normalized response')
        end
    end
    
    %% Validate functional topography
    figure
    
    % select random units
    nExamps = 3;
    unitInds = randsample(size(Rinit,4),nExamps);
    
    % Find K nearest neighbors in functional space
    neighbors = knnsearch(Y,Y(unitInds,:),'K',Kneighbors);
    
    % For each neuron, plot PSTH of nearest neighbors
    for exUnit = 1:nExamps
        subplot(nExamps,2,1+(exUnit-1)*2)
        for ni = 1:Kneighbors
            if dynCohCollate
                zR = (squeeze(Rdyn(:,5,neighbors(exUnit,ni)))-mean(squeeze(Rdyn(:,5,neighbors(exUnit,ni)))))/std(squeeze(Rdyn(:,5,neighbors(exUnit,ni))));
                plot(dynCoh.neuron_t(dynCoh.neuron_t<=1350),zR,'Color',colors(ni,:))
                hold on
            elseif initCohCollate
                zR = (squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))-mean(squeeze(Rinit(:,3,3,neighbors(exUnit,ni)))))/std(squeeze(Rinit(:,3,3,neighbors(exUnit,ni))));
                plot(initCoh.neuron_t(initCoh.neuron_t<=1350),zR(initCoh.neuron_t<=1350),'Color',colors(ni,:))
                hold on
            end
        end
        xlabel('Time from motion onset (ms)')
        ylabel('z-score')
        
        subplot(nExamps,2,2+(exUnit-1)*2)
        plot(Y(:,1),Y(:,2),'k.')
        hold on
        for ni = 1:Kneighbors
            plot(Y(neighbors(exUnit,ni),1),Y(neighbors(exUnit,ni),2),'o','Color',colors(ni,:))
        end
    end
    
    %% topography
    hChamberTopology = figure;
    for subjecti = 1:length(subjects)
        subplot(1,length(subjects),subjecti)
        randScale = 0.08;
        locations2 = locations(cellID(:,1,4) == subjecti,:);
        locations2(locations2(:,1)>1,:) = [locations2(locations2(:,1)>1,2), locations2(locations2(:,1)>1,1), locations2(locations2(:,1)>1,3)];     % correction for mapping error between tetrode and vprobe recording set ups
        locationsRand = locations2 + ...
            [randScale*nanstd(locations2(:,1))*randn(size(locations2,1),1), ...
            randScale*nanstd(locations2(:,2))*randn(size(locations2,1),1), ...
            0*nanstd(locations2(:,3))*randn(size(locations2,1),1)];                  % Add randomness to a-p and m-l locations to make
        for uniti = 1:size(locations2,1)
            plot3(locationsRand(uniti,1),locationsRand(uniti,2),locationsRand(uniti,3)/1000,...
                'o','Color',[0 0 0],'MarkerSize',6,'MarkerFaceColor',colorWheel(idx(uniti),:,:));
            hold on
        end
        grid on
        axis equal
        xlabel('Anterior/postieror (mm)')
        ylabel('Medial/lateral (mm)')
        zlabel('Depth (mm)')
        title([subjects{subjecti} ' recoding locations'])
    end
    
    %% Functional vs physical topography
    figure('Name','Functional vs physical topography','Position',[1153 924 1373 405])
    for subjecti = 1:length(subjects)
        subplot(length(subjects),2,1+2*(subjecti-1))
        samps = randsample(length(euclidLoc{subjecti}),500,false);
        plot(euclidLoc{subjecti}(samps),euclidY{subjecti}(samps),subjectMarkers{subjecti})
        hold on
        ax = axis;
        text((ax(2)-ax(1))*0.8+ax(1),(ax(4)-ax(3))*0.9+ax(3),['R = ' num2str(locCor(1,2))])
        text((ax(2)-ax(1))*0.8+ax(1),(ax(4)-ax(3))*0.85+ax(3),['p = ' num2str(locCorP(1,2))])
        xlabel('Physical distance (mm)')
        ylabel('Functional distance (a.u.)')
        title(['Subject ' subjects{subjecti}])
        
        subplot(length(subjects),2,2+2*(subjecti-1))
        for ni = 1:size(nnIdx{subjecti},1)
            mdist(ni) = mean(pdist(Y(nnIdx{subjecti}(ni,:),:)));
            randIdx = [ni; randsample(size(Y,1),Kneighbors-1)];
            mdistRand(ni) = mean(pdist(Y(randIdx,:)));
        end
        histogram(mdist,linspace(min([mdist,mdistRand]),max([mdist,mdistRand]),50))
        hold on
        histogram(mdistRand,linspace(min([mdist,mdistRand]),max([mdist,mdistRand]),50))
        xlabel('Functional distance')
        ylabel('N')
        legend({[num2str(Kneighbors) ' nearest in chamber'],'Random'})
        title(['Subject ' subjects{subjecti}])
    end
    
    %% Heat map of population firing rates, sorted by embedding location 
    if initCohCollate
        hInitCohRatesHeatMap = figure('Name','Firing rate heat map, initCoh','Position',[864 67 560 1169]);
        cax = [Inf,-Inf];
        RinitNorm = Rinit./max(Rinit,[],[1,2,3]);
        for speedi = 1:length(speeds)
            for cohi = 1:length(cohs)
                subplot(length(speeds),length(cohs),cohi + (speedi-1)*length(cohs))
                Rtemp = squeeze(RinitNorm(:,speedi,cohi,:));
                %             Rtemp = Rtemp./max(Rtemp,[],1);
                imagesc(initCoh.neuron_t,1:size(Rtemp,2),Rtemp(:,thetaSort)')
                xlabel('Time from motion onset (ms)')
                ylabel('Neuron #')
                
                caxTemp = caxis;
                cax(1) = min([caxTemp(1) cax(1)]);
                cax(2) = max([caxTemp(2) cax(2)]);
            end
        end
        for speedi = 1:length(speeds)
            for cohi = 1:length(cohs)
                subplot(length(speeds),length(cohs),cohi + (speedi-1)*length(cohs))
                caxis(cax)
                set(gca,'TickDir','out')
            end
        end
    end
    
    if dynCohCollate
        hDynCohRatesHeatMap = figure('Name','Firing rate heat map, dynCoh','Position',[864 67 560*5/3 1169/3]);
        cax = [Inf,-Inf];
        RdynNorm = Rdyn./max(Rdyn,[],[1,2]);
        for seqi = 1:length(sequences)
            subplot(1,length(sequences),seqi)
            Rtemp = squeeze(RdynNorm(:,seqi,:));
            %             Rtemp = Rtemp./max(Rtemp,[],1);
            imagesc(dynCoh.neuron_t(dynCoh.neuron_t<=1350),1:size(Rtemp,2),Rtemp(:,thetaSort)')
            xlabel('Time from motion onset (ms)')
            ylabel('Neuron #')
            
            caxTemp = caxis;
            cax(1) = min([caxTemp(1) cax(1)]);
            cax(2) = max([caxTemp(2) cax(2)]);
        end
        for seqi = 1:length(sequences)
            subplot(1,length(sequences),seqi)
            caxis(cax)
            set(gca,'TickDir','out')
        end
    end
    
    
    %% Save figures
    if saveFigures
        
        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' [subjects{:}] ...
            '/neuronTyping/' datestr(now,'yyyymmdd')];
        if ~exist(saveLocation,'dir')
            mkdir(saveLocation)
        end
        
        for uniti = 1:size(exUnits,1)
            exIndex = exUnits(uniti,:);
            cellID2 = squeeze(cellID(:,1,1:2));
            listIndex = find(ismember(cellID2, exIndex, 'rows'));
            savefig(hExUnits(uniti),[saveLocation '/exNeuron_' num2str(exIndex(1)) '_' num2str(exIndex(2)) '_' subjects{cellID(listIndex,1,4)} '.fig'])
        end
        for i = 1:NumClusters
            savefig(gvrh(i), [saveLocation '/gainVsClusterFiringRates' num2str(i) '.fig'])
            savefig(hClusterPSTH(i), [saveLocation '/clusterPSTH' num2str(i) '.fig'])
        end
        savefig(hGainDecodingFromClusters , [saveLocation '/gainDecodingClusters.fig'])
        savefig(hFunctionalTopology ,[saveLocation '/functionalTopology.fig'])
        savefig(hGainTopology,[saveLocation '/gainTopology.fig'])
        savefig(hChamberTopology, [saveLocation '/chamberTopology.fig'])
        if initCohCollate
            savefig(hClusterPSTHsInitCoh ,[saveLocation '/clusterPSTHsInitCoh.fig'])
            savefig(hInitCohRatesHeatMap ,[saveLocation '/initCohRatesHeatMap.fig'])
        end
        if dynCohCollate
            savefig(hDynCohRatesHeatMap ,[saveLocation '/dynCohRatesHeatMap.fig'])  
            savefig(hClusterPSTHsDynCoh ,[saveLocation '/clusterPSTHsDynCoh.fig'])
        end
        
    end