function Beta = eggerMTdecoder(varargin)
%% MTtoFEFregression
%
%
%%

%% Defaults
plotOpts_default.On = false;
speedPrefOpts_default.tWin = [40,120];
speedPrefOpts_default.P0 = [16,1];
speedPrefOpts_default.ub = [254,128];
speedPrefOpts_default.lb = [0, 0];
speedPrefOpts_default.c = NaN;
speedPrefOpts_default.s = NaN;
speedPrefOpts_default.d = 0;
compToBehavior_default.On = false;
theoretical_default.weightTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';
simulateMT_default.On = false;
estimatorPathwayNonlinearity_default = @(x)(x);

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'epsilon',0)
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'estimatorPathwayNonlinearity',estimatorPathwayNonlinearity_default)
addParameter(Parser,'tWin',[0 900])
addParameter(Parser,'rankN',80)
addParameter(Parser,'ridgeLambda',logspace(-1,12,10))
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'compToBehavior',compToBehavior_default)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'simulateMT',simulateMT_default)
addParameter(Parser,'normalizeMT',true)
addParameter(Parser,'plotOpts',plotOpts_default)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,varargin{:})

objectFile = Parser.Results.objectFile;
epsilon = Parser.Results.epsilon;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
estimatorPathwayNonlinearity = Parser.Results.estimatorPathwayNonlinearity;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
compToBehavior = Parser.Results.compToBehavior;
theoretical = Parser.Results.theoretical;
simulateMT = Parser.Results.simulateMT;
normalizeMT = Parser.Results.normalizeMT;
plotOpts = Parser.Results.plotOpts;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;


%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit
if simulateMT.On
    mtNeuron_t = -500:900;
    speeds = speedsFEF;
    [MT, spref, cohs] = simulateMTdata(mtResults,mtNeuron_t,simulateMT.modelN,...
        'speedsMT',speeds,'removeBaseline',simulateMT.removeBaseline,'gaussianApprox',simulateMT.gaussianApprox);
    swidth = 1.2*ones(size(spref));
else
    [MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
        'speedsMT',speeds,'cohsMT',cohs,'directionsMT',directionsMT,...
        'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
        'speedPrefOpts',speedPrefOpts);
    
    
    mtNeuron_t = mtResults.mt{1}.neuron_t;
end

% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);

mt = permute(interpolatedR,[4,1,2,3]);
mt = mt(prod(any(isnan(mt),2),[3,4])==0,:,:,:);

t = mtNeuron_t;
mt_z = (mt-mean(mt(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(mt(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);


sp = spref(~isnan(interpolatedR(1,1,1,:)));
sw = swidth(~isnan(interpolatedR(1,1,1,:)));
[~,sortInd] = sort(sp);
[n,x] = histcounts(sp,10,'Normalization','pdf');
p = interp1(x(2:end)-(x(2)-x(1))/2,n,sp,'linear','extrap');
ip = mean(p)./p;

%% decoding Weights
switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(sp)') ones(size(sp'))];
        
    case 'optimal'
        % More complicated: 'optimal' decoder assuming zero correlations:
        % df/ds*I/(df/ds'*I*df/ds)
        s0 = speedsFEF(speedsFEF==10);
        sprefTemp = sp;
        %     sprefTemp(sprefTemp < 1) = 1;
        swidthTemp = sw;
        %     swidthTemp(swidthTemp > 10) = 10;
        df = -log2(s0./sprefTemp)./(s0.*sw*log(2)).*exp(log2(s0./sprefTemp).^2./(2*sw.^2));
        uOpt = df*inv(eye(length(sprefTemp)))/(df*inv(eye(length(sprefTemp)))*df');
        uOpt(uOpt<-0.05) = min(uOpt(uOpt>-0.05));
        uOptNorm = uOpt/norm(uOpt);
        Atheory = [uOpt'-min(uOpt), ones(size(uOpt'))];
end

% Normalize 
%Atheory = Atheory./vecnorm(Atheory);

%% Calculate estimate and gain pathway outputs
if normalizeMT
    mtApplied = mt_z;
else
    mtApplied = mt;
end
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        shat(:,:,si,ci) = Atheory(:,:)'*mtApplied(:,:,si,ci) ./ (epsilon + sum(mtApplied(:,:,si,ci),1)); %log2(speedsFEF(si))*ones(size(Atheory'*mt(:,:,si,ci))); %
        gain(:,:,si,ci) = Atheory(:,:)'*mtApplied(:,:,si,ci);
    end
end

shat = estimatorPathwayNonlinearity(shat);
ehat = gain.*(shat);
c = median(speedsFEF)/ehat(1,t==80,speedsFEF==median(speedsFEF),cohsFEF==100);

gain = gain*c;
ehat = gain.*(shat);

%% Saving
if saveResults
    saveLocation = ['/home/seth/Projects/DynamicCoherencePhysiology/MTmodel'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/EggerDecoder' datestr(now,'yyyymmdd')],'-v7.3')
end
    

%% Plotting
if saveFigures
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/MTmodel' ...
        '/EggerDecoder/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
end

%% Plot results
if plotOpts.On
    
    %% Decoder output over time
    hOutputOverTime = figure;
    for ci = 1:length(cohsFEF)
        subplot(3,length(cohsFEF),ci)
        for si = 1:length(cohsFEF)
            plot(t,shat(1,:,si,ci),'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Estimate pathway (deg/s)')
        set(gca,'TickDir','out')
        ax(ci,:,1) = axis;
        
        subplot(3,length(cohsFEF),3+ci)
        for si = 1:length(cohsFEF)
            plot(t,gain(1,:,si,ci),'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Gain')
        set(gca,'TickDir','out')
        ax(ci,:,2) = axis;
        
        subplot(3,length(cohsFEF),6+ci)
        for si = 1:length(cohsFEF)
            plot(t,ehat(1,:,si,ci),'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel('Eye speed (deg/s)')
        set(gca,'TickDir','out')
        ax(ci,:,3) = axis;
        
    end
    
    for ci = 1:length(cohsFEF)
        subplot(3,length(cohsFEF),ci)
        if max(ax(:,4,1)) < 3*max(speedsFEF)
            axis([min(ax(:,1,1)) max(ax(:,2,1)) min(ax(:,3,1)) max(ax(:,4,1))])
        else
            axis([min(ax(:,1,1)) max(ax(:,2,1)) min(ax(:,3,1)) 3*max(speedsFEF)])
        end
        
        subplot(3,length(cohsFEF),3+ci)
        axis([min(ax(:,1,2)) max(ax(:,2,2)) min(ax(:,3,2)) max(ax(:,4,2))])
        
        subplot(3,length(cohsFEF),6+ci)
        axis([min(ax(:,1,3)) max(ax(:,2,3)) min(ax(:,3,3)) max(ax(:,4,3))])
    end
    
    if saveFigures
        savefig(hOutputOverTime,[saveLocation '/decoderOutputOverTime.fig'])
    end
    
    %% MT response
    hMTInitiation = figure('Name','Population response, intiation');
    lims = [Inf -Inf];
    edges = logspace(0,log10(128),50);
    [~,bin] = histc(sp,logspace(0,log10(128),100));
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            tempMT = mtApplied(:,t==80,si,ci);
            tempMTave = nan(length(edges),1);
            for bi = 1:length(edges)
                tempMTave(bi) = mean(tempMT(bin == bi));
            end
            subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
%             semilogx(sp(:),tempMT,'o','Color',speedColors(si,:))
            s = scatter(sp(:),tempMT,80,repmat(speedColors(si,:),[length(sp),1]),'filled');
            s.MarkerFaceAlpha = 0.1;
            hold on
%             semilogx(edges + (edges(2)-edges(1))/2,smooth(tempMTave),...
%                 'Color',[0 0 0],'LineWidth',2)
            xlabel('Preferred speed (deg/s)')
            ylabel('Pop response')
            axis tight
            limTemp = ylim;
            lims(1) = min([lims(1),limTemp(1)]);
            lims(2) = max([lims(2),limTemp(2)]);
        end
    end
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
            set(gca,'XScale','log','TickDir','out')
            ylim(lims)
            hold on
            plotVertical(shat(1,t==80,si,ci));
        end
    end
    
    if saveFigures
        savefig(hMTInitiation,[saveLocation '/MTpop_Initiation.fig'])
    end
        
    hMTSteadyState = figure('Name','Population response, steady state tracking');
    lims = [Inf -Inf];
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
%             semilogx(sp(:),mtApplied(:,t==600,si,ci),'o','Color',speedColors(si,:))
            s2 = scatter(sp(:),mtApplied(:,t==600,si,ci),80,repmat(speedColors(si,:),[length(sp),1]),'filled');
            s2.MarkerFaceAlpha = 0.1;
            xlabel('Preferred speed (deg/s)')
            ylabel('Pop response')
            axis tight
            limTemp = ylim;
            lims(1) = min([lims(1),limTemp(1)]);
            lims(2) = max([lims(2),limTemp(2)]);
        end
    end
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(speedsFEF),length(cohsFEF),ci+(si-1)*length(cohsFEF))
            set(gca,'XScale','log','TickDir','out')
            ylim(lims)
            hold on
            plotVertical(shat(1,t==600,si,ci));
        end
    end
    
    if saveFigures
        savefig(hMTSteadyState,[saveLocation '/MTpop_SteadyState.fig'])
    end
    
    
    %% Outputs during initiaton
    hInitiationOutput = figure('Name','Pathway outputs','Position',[1923 15 560 1286]);
    subplot(4,1,1)
    for ci = 1:length(cohsFEF)
        plot(speedsFEF,2.^squeeze(shat(1,t==80,:,ci)),'o-',...
            'Color',ones(1,3)-cohsFEF(ci)/100,'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100,...
            'DisplayName',[num2str(cohsFEF(ci)) '%'])
        hold on
    end
    axis square
    axis([4 22 4 22])
    plotUnity;
    xlabel('Target speed (deg/s)')
    ylabel('2^{(Esimate pathway)} (deg/s)')
    set(gca,'TickDir','out')
    
    subplot(4,1,2)
    for ci = 1:length(cohsFEF)
        plot(speedsFEF,squeeze(shat(1,t==80,:,ci)),'o-',...
            'Color',ones(1,3)-cohsFEF(ci)/100,'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100,...
            'DisplayName',[num2str(cohsFEF(ci)) '%'])
        hold on
    end
    axis square
    xlim([4 22])
    plotUnity;
    xlabel('Target speed (deg/s)')
    ylabel('Estimate pathway')
    set(gca,'TickDir','out')
    
    subplot(4,1,3)
    for ci = 1:length(cohsFEF)
        plot(speedsFEF,squeeze(gain(1,t==80,:,ci)),'o-',...
            'Color',ones(1,3)-cohsFEF(ci)/100,'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100,...
            'DisplayName',[num2str(cohsFEF(ci)) '%'])
        hold on
    end
    axis square
    xlim([4 22])
    xlabel('Target speed (deg/s)')
    ylabel('FEF pathway')
    set(gca,'TickDir','out')
    
    subplot(4,1,4)
    for ci = 1:length(cohsFEF)
        plot(speedsFEF,squeeze(ehat(1,t==80,:,ci)),'o-',...
            'Color',ones(1,3)-cohsFEF(ci)/100,'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100,...
            'DisplayName',[num2str(cohsFEF(ci)) '%'])
        hold on
    end
    axis square
    axis([4 22 4 22])
    plotUnity;
    xlabel('Target speed (deg/s)')
    ylabel('Estimated eye speed (deg/s)')
    set(gca,'TickDir','out')
    
    if saveFigures
        savefig(hInitiationOutput,[saveLocation '/decoderInitationOutput.fig'])
    end
    
    %% Comparison to behavior
    if compToBehavior.On
        hSpeeds = figure;
        hGains = figure;
        for filei = 1:length(compToBehavior.file)
            temp = load(compToBehavior.file{filei},'init');
            eye_t = temp.init.t;
            eyeSpeed = temp.init.eye.mean;
            slips = -squeeze(temp.init.eye.mean(temp.init.t == 750,:,:)) + speedsFEF';
            initGain = (temp.init.eye.pert.res - temp.init.eye.pert.resControl)./(0.4*repmat(speedsFEF',[1,3,3])); 
            clear temp
            if compToBehavior.applyLogrithmicEstimatorCorrection
                gBehavior(:,:,1) = initGain(:,:,2).*speedsFEF'*0.4./log2(1.4);
                gBehavior(:,:,2) = initGain(:,:,3).*speedsFEF'*0.4./log2(1+0.4.*speedsFEF'./slips);
            else
                gBehavior = initGain(:,:,2:3);
            end
            ehatBehavioral = squeeze(shat(1,t==80,:,:)).*gBehavior;
            ehatBehavioral_z = (ehatBehavioral - repmat(mean(ehatBehavioral,[1,2]),[length(speedsFEF),length(cohsFEF),1]))./ ...
                repmat(std(ehatBehavioral,[],[1,2]),[length(speedsFEF),length(cohsFEF),1]);
            
            figure(hSpeeds)
            subplot(length(compToBehavior.file),3,1+(filei-1)*3)
            for ci = 1:length(cohsFEF)
                eyeTemp = (eyeSpeed(eye_t==150,:,ci)-...
                    mean(eyeSpeed(eye_t==150,:,:),[2,3]))./...
                    std(eyeSpeed(eye_t==150,:,:),[],[2,3]);
                plot(speedsFEF,eyeTemp,'o-',...
                    'Color',ones(1,3)-cohsFEF(ci)/100)
                hold on
                plot(speedsFEF,(squeeze(ehat(1,t==80,:,ci))-mean(ehat(1,t==80,:,:),[3,4]))./std(ehat(1,t==80,:,:),[],[3,4]),...
                    'x','Color',ones(1,3)-cohsFEF(ci)/100)
%                 plot(speedsFEF,ehatBehavioral_z(:,ci,1),...
%                     'd','Color',ones(1,3)-cohsFEF(ci)/100,'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100)
            end
            xlabel('Target speed (deg/s)')
            ylabel('z-scored eye speed')
            set(gca,'TickDir','out')
            
            subplot(length(compToBehavior.file),3,2+(filei-1)*3)
            for ci = 1:length(cohsFEF)
                eyeTemp = (eyeSpeed(eye_t==150,:,ci)-...
                    mean(eyeSpeed(eye_t==150,:,:),[2,3]))./...
                    std(eyeSpeed(eye_t==150,:,:),[],[2,3]);
                plot(eyeTemp,...
                    (squeeze(ehat(1,t==80,:,ci))-mean(ehat(1,t==80,:,:),[3,4]))./std(ehat(1,t==80,:,:),[],[3,4]),...
                    'x-','Color',ones(1,3)-cohsFEF(ci)/100)
                hold on
%                 plot(eyeTemp,...
%                     ehatBehavioral_z(:,ci,1),...
%                     'd-','Color',ones(1,3)-cohsFEF(ci)/100)
            end
            plotUnity();
            axis square
            xlabel('z-scored eye speed')
            ylabel('z-scored predictions')
            set(gca,'TickDir','out')
            
            subplot(length(compToBehavior.file),3,3+(filei-1)*3)
            for ci = 1:length(cohsFEF)
                eyeTemp = (eyeSpeed(eye_t==750,:,ci)-...
                    mean(eyeSpeed(eye_t==150,:,:),[2,3]))./...
                    std(eyeSpeed(eye_t==150,:,:),[],[2,3]);
                plot(eyeTemp,...
                    (squeeze(ehat(1,t==750,:,ci))-mean(ehat(1,t==80,:,:),[3,4]))./std(ehat(1,t==80,:,:),[],[3,4]),...
                    'x-','Color',ones(1,3)-cohsFEF(ci)/100)
                hold on
%                 plot(eyeTemp,...
%                     ehatBehavioral_z(:,ci,1),...
%                     'd-','Color',ones(1,3)-cohsFEF(ci)/100)
            end
            plotUnity();
            axis square
            xlabel('z-scored eye speed')
            ylabel('z-scored predictions')
            set(gca,'TickDir','out')
            
            figure(hGains)
            subplot(2,2,filei)
            for ci = 1:length(cohsFEF)
                plot(gBehavior(:,ci,1),...
                    squeeze(gain(1,t==80,:,ci)),...
                    'o-','Color',ones(1,3)-cohsFEF(ci)/100,...
                    'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100)
                hold on
            end
            xlabel('Measured gain')
            ylabel('Predicted gain')
            set(gca,'TickDir','out')
            ax = axis;
            for ci = 1:length(cohsFEF)
                gtemp = gain(1,t==80,:,ci);
                rtemp = corrcoef(gBehavior(:,ci,1),gtemp(:));
                text(0.05*(ax(2)-ax(1))+ax(1),(0.95-0.1*(ci-1))*(ax(4)-ax(3))+ax(3),['r = ' num2str(rtemp(1,2))],'Color',ones(1,3)-cohsFEF(ci)/100)
            end
            
            subplot(2,2,filei+2)
            for ci = 1:length(cohsFEF)
                plot(gBehavior(:,ci,2),...
                    squeeze(gain(1,t==750,:,ci)),...
                    'o-','Color',ones(1,3)-cohsFEF(ci)/100,...
                    'MarkerFaceColor',ones(1,3)-cohsFEF(ci)/100)
                hold on
            end
            xlabel('Measured gain')
            ylabel('Predicted gain')
            set(gca,'TickDir','out')
            ax = axis;
            for ci = 1:length(cohsFEF)
                gtemp = gain(1,t==750,:,ci);
                rtemp = corrcoef(gBehavior(:,ci,2),gtemp(:));
                text(0.05*(ax(2)-ax(1))+ax(1),(0.95-0.1*(ci-1))*(ax(4)-ax(3))+ax(3),['r = ' num2str(rtemp(1,2))],'Color',ones(1,3)-cohsFEF(ci)/100)
            end
            
        end
        
        
        if saveFigures
            savefig(hSpeeds,[saveLocation '/decoderVsSpeeds.fig'])
        end
        
        if saveFigures
            savefig(hGains,[saveLocation '/decoderVsGains.fig'])
        end
        
    end
end
