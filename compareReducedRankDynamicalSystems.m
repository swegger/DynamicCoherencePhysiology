function compareReducedRankDynamicalSystems(subject,varargin)
%%
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'fitTypes',{'fitSigmas','fitSigmasWithDummy','fitSigmasExtraInput'})
addParameter(Parser,'generationDate',[])
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20,60,100])
addParameter(Parser,'saveResults',false)
addParameter(Parser,'saveFigures',false)

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
fitTypes = Parser.Results.fitTypes;
generationDate = Parser.Results.generationDate;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
saveResults = Parser.Results.saveResults;
saveFigures = Parser.Results.saveFigures;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Load results
model_cc.eyeSpeed = nan(length(fitTypes),300);
model_cc.gain = nan(length(fitTypes),300);
model_mse = nan([1001,300,2,length(fitTypes)]);
for fitTypei = 1:length(fitTypes)
    
    disp(['Fit type ' fitTypes{fitTypei}])
    
    sourcePath = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' subject '/fitReducedRankModel/' fitTypes{fitTypei} '/'];
    if ~isempty(generationDate)
        sourcePattern = ['*_' generationDate '.mat'];
    else
        sourcePattern = [];
    end
    files = dir([sourcePath, sourcePattern]);
        
    for filei = 1:length(files)
        try
            temp = load([sourcePath, files(filei).name],'N','speeds','meanEyeSpeed','eye_t','Rhat','t','initGain','RtoFit','dimNames','modelFEF');
            num = regexp(files(filei).name,'_\d+_','match');
            num = str2num(num{1}(num{1}~='_'));
            tempData = [];
            for si = 1:length(temp.speeds)
                eyeSpeedTemp = squeeze(nanmean(temp.meanEyeSpeed(temp.eye_t >= 750 & temp.eye_t <= 750,:,:),1));
                targetedDim = find(strcmp(temp.dimNames,'Speed'));
                initRatesTempSpeed = squeeze(nanmean(temp.Rhat(temp.t >= 700 & temp.t <= 800,si,:,targetedDim),1));
                targetedDim = find(strcmp(temp.dimNames,'Gain'));
                initRatesTempGain = squeeze(nanmean(temp.Rhat(temp.t >= 700 & temp.t <= 800,si,:,targetedDim),1));
                tempData = [tempData; eyeSpeedTemp(si,:)', temp.initGain(si,:,3)', initRatesTempSpeed(:) initRatesTempGain(:)];
            end
            
            for si = 1:length(temp.speeds)                
                eyeSpeedTemp = squeeze(nanmean(temp.meanEyeSpeed(temp.eye_t >= 150 & temp.eye_t <= 150,:,:),1));
                targetedDim = find(strcmp(temp.dimNames,'Speed'));
                initRatesTempSpeed = squeeze(nanmean(temp.Rhat(temp.t >= 100 & temp.t <= 200,si,:,targetedDim),1));
                targetedDim = find(strcmp(temp.dimNames,'Gain'));
                initRatesTempGain = squeeze(nanmean(temp.Rhat(temp.t >= 100 & temp.t <= 200,si,:,targetedDim),1));
                tempData = [tempData; eyeSpeedTemp(si,:)', temp.initGain(si,:,2)', initRatesTempSpeed(:) initRatesTempGain(:)];
            end
            rvals = corrcoef(tempData);
            model_cc.eyeSpeed(fitTypei,num) = rvals(1,3);
            model_cc.gain(fitTypei,num) = rvals(2,4);
            
            model_mse(:,num,:,fitTypei) = mean((temp.RtoFit-temp.Rhat(:,:,:,length(temp.dimNames))).^2,[2,3]);
            
        catch
            disp([sourcePath, files(filei).name])
        end
            

    end
    disp(num2str(sum(~isnan(model_cc.gain(fitTypei,:)))))
end

%% Saving
if saveResults
        
        saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/' subject ...
            '/ReducedRankModel'];
        if ~exist(saveLocation,'dir')
            mkdir(saveLocation)
        end
            save([saveLocation '/modelComparison' datestr(now,'yyyymmdd')],'-v7.3')
    
end


%% Plot results


%% Summary of behavioral predictive ability
hsummaryHist = figure('Name','Histogram of model output, behavioral data predictive power');
subplot(1,2,1)
for fitTypei = 1:length(fitTypes)
    histogram(model_cc.eyeSpeed(fitTypei,:),linspace(-1,1,21),...
        'DisplayName',[fitTypes{fitTypei}]);
    hold on
end
xlabel('Model and eye speed R')
ylabel('# of fits')

subplot(1,2,2)
for fitTypei = 1:length(fitTypes)
    histogram(model_cc.gain(fitTypei,:),linspace(-1,1,21),...
        'DisplayName',[fitTypes{fitTypei}]);
    hold on
end
xlabel('Model and gain R')
ylabel('# of fits')

hspeedVgainPerformance = figure('Name','Speed vs gain represenation for same MT inputs');
for fitTypei = 1:length(fitTypes)
    subplot(1,length(fitTypes),fitTypei)
    plot(model_cc.eyeSpeed(fitTypei,:),model_cc.gain(fitTypei,:),'o')
    xlabel('Speed prediction performance')
    ylabel('Gain prediction performance')
    axis([-1 1 -1 1])
    axis square
end
    
hmodelVmodelPerformance = figure('Name','Across model correlatoins for same MT inputs');
ind = 0;
for i = 1:length(fitTypes)
    for j = 1:length(fitTypes)
        if j > i
            ind = ind+1;
            subplot(1,length(fitTypes)*(length(fitTypes)-1)/2,ind)
            plot(model_cc.eyeSpeed(i,:),model_cc.eyeSpeed(j,:),'o',...
                'DisplayName','Eye speed')
            hold on
            plot(model_cc.gain(i,:),model_cc.gain(j,:),'o',...
                'DisplayName','Gain')
            xlabel(['R values ' fitTypes{i}])
            ylabel(['R values ' fitTypes{j}])
            axis([-1 1 -1 1])
            axis square
            plotUnity;
        end
    end
end


%% Summary of target function prediction
figure('Name','RMSE over time, across conditions')
for fitTypei = 1:length(fitTypes)
    for di = 1:length(temp.dimNames)
        acceptVec = any(isnan(model_mse(:,:,di,fitTypei)));
        subplot(length(temp.dimNames),length(fitTypes),fitTypei+(di-1)*length(fitTypes))
        imagesc(1:sum(~acceptVec),linspace(-100,900,1001),sqrt(model_mse(:,~acceptVec,di,fitTypei)))
        xlabel('MT population simulation')
        ylabel('Time from motion onset (ms)')
        caxis([min(sqrt(model_mse(:,:,di,:)),[],[1,2,4]) max(sqrt(model_mse(:,:,di,:)),[],[1,2,4])])
    end
end
    
%% Example model behavior for most predictive models of each type

for fitTypei = 1:length(fitTypes)
    
    % Load data
    sourcePath = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' subject '/fitReducedRankModel/' fitTypes{fitTypei} '/'];
    [~,num] = max(model_cc.eyeSpeed(fitTypei,:).^2 + model_cc.eyeSpeed(fitTypei,:).^2);
    nums(fitTypei) = num;
    filename = ['fitReducedRankModel_' fitTypes{fitTypei} '_' num2str(num) '_' generationDate '.mat'];
    temp = load([sourcePath, filename],'N','speeds','meanEyeSpeed','eye_t',...
        'Rhat','t','initGain','RtoFit','dimNames','modelFEF','theoreticalInput');
    temp.M = size(temp.theoreticalInput,1);
    
    hModelandData(fitTypei) = figure('Name',['Predicted response along targeted dimensions, ' fitTypes{fitTypei}],'Position',[63 169 1606 1079]);
    for dimi = 1:temp.M
        subplot(temp.N+temp.M,2,2*(dimi-1)+1)
        for speedi = 1:length(speedsFEF)
            for cohi = 1:length(cohsFEF)
                plot(temp.t,temp.theoreticalInput(dimi,:,speedi,cohi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100],...
                    'DisplayName',['Speed = ' num2str(speedsFEF(speedi)) ', Coh = ' num2str(cohsFEF(cohi)), ' model = ' fitTypes{fitTypei}])
                hold on
            end
        end
        axis tight
        title('Input from MT')
        xlabel('Time from motion onset (ms)')
        ylabel(['Input ' num2str(dimi)])
    end
    
    
    for dimi = 1:temp.N
        subplot(temp.N+temp.M,2,2*(dimi-1)+2*temp.M+1)
        for speedi = 1:length(speedsFEF)
            for cohi = 1:length(cohsFEF)
                plot(temp.t,temp.Rhat(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100],...
                    'DisplayName',['Speed = ' num2str(speedsFEF(speedi)) ', Coh = ' num2str(cohsFEF(cohi)), ' model = ' fitTypes{fitTypei}])
                hold on
            end
        end
        axis tight
        title('Fit predictions')
        xlabel('Time from motion onset (ms)')
        if dimi > length(temp.dimNames)
            ylabel(['Dummy related activity (a.u.)'])
        else
            ylabel([temp.dimNames{dimi} ' related activity (a.u.)'])
        end
        
        subplot(temp.N+temp.M,2,2*(dimi-1)+2*temp.M+2)
        if dimi <= length(temp.dimNames)
            for speedi = 1:length(speedsFEF)
                for cohi = 1:length(cohsFEF)
                    plot(temp.t,temp.RtoFit(:,speedi,cohi,dimi),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100],...
                        'DisplayName',['Speed = ' num2str(speedsFEF(speedi)) ', Coh = ' num2str(cohsFEF(cohi)), ' model = ' fitTypes{fitTypei}])
                    hold on
                end
            end
            axis tight
            title('Data')
            xlabel('Time from motion onset (ms)')
            ylabel([temp.dimNames{dimi} ' related activity (a.u.)'])
        end
    end
    
    compLims = repmat([Inf,-Inf],[temp.N,1]);
    for dimi = 1:temp.N
        for compi = 1:2
            subplot(temp.N+temp.M,2,2*(dimi-1)+2*temp.M+compi)
            limTemp = ylim;
            compLims(dimi,1) = min([compLims(dimi,1),limTemp(1)]);
            compLims(dimi,2) = max([compLims(dimi,2),limTemp(2)]);
        end
    end
    for dimi = 1:temp.N
        for compi = 1:2
            subplot(temp.N+temp.M,2,2*(dimi-1)+2*temp.M+compi)
            ylim(compLims(dimi,:))
        end
    end
        
    gvth(fitTypei) = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh), ' fitTypes{fitTypei}],'Position',[1456 59 1070 1263]);
    for targetedDim = 1:temp.N
        tempData = [];
        subplot(temp.N,2,2*targetedDim-1)
        for si = 1:length(speedsFEF)
            eyeSpeedTemp = squeeze(nanmean(temp.meanEyeSpeed(temp.eye_t >= 750 & temp.eye_t <= 750,:,:),1));
            initRatesTemp = squeeze(nanmean(temp.Rhat(temp.t >= 700 & temp.t <= 800,si,:,targetedDim),1));
            plot(eyeSpeedTemp(si,:),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si))])
            hold on
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        for si = 1:length(speedsFEF)
            eyeSpeedTemp = squeeze(nanmean(temp.meanEyeSpeed(temp.eye_t >= 150 & temp.eye_t <= 150,:,:),1));
            initRatesTemp = squeeze(nanmean(temp.Rhat(temp.t >= 100 & temp.t <= 200,si,:,targetedDim),1));
            plot(eyeSpeedTemp(si,:),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Speed = ' num2str(speedsFEF(si))])
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        
        speed_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
        axis square
        ax = axis;
        text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(speed_projection_cc(1,2).^2)])
        xlabel('Eye speed (deg/s)')
        if targetedDim <= length(temp.dimNames)
            ylabel(['Projection on ' temp.dimNames{targetedDim} ' dimension'])
        end
        
        tempData = [];
        subplot(temp.N,2,2*targetedDim)
        for si = 1:length(speedsFEF)
            initRatesTemp = squeeze(nanmean(temp.Rhat(temp.t >= 700 & temp.t <= 800,si,:,targetedDim),1));
            plot(squeeze(temp.initGain(si,:,3)),initRatesTemp,...
                'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Steady State pursuit, Speed = ' num2str(speedsFEF(speedi))])
            hold on
            tempData = [tempData; temp.initGain(si,:,3)',initRatesTemp(:)];
            initRatesTemp = squeeze(nanmean(temp.Rhat(temp.t >= 100 & temp.t <= 200,si,:,targetedDim),1));
            plot(squeeze(temp.initGain(si,:,2)),initRatesTemp,...
                'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
                'DisplayName',['Pursuit initiation, Speed = ' num2str(speedsFEF(speedi))])
            tempData = [tempData; temp.initGain(si,:,2)',initRatesTemp(:)];
        end
        gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
        axis square
        ax = axis;
        text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
        xlabel('Behavioral gain (unitless)')
        if targetedDim <= length(temp.dimNames)
            ylabel(['Projection on ' temp.dimNames{targetedDim} ' dimension'])
        end
    end
    
end

%% Save figures
if saveFigures
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' subject ...
        '/fitReducedRankModel/modelComparison/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    
    savefig(hsummaryHist,[saveLocation '/performanceSummary.fig'])
    savefig(hmodelVmodelPerformance,[saveLocation '/pairedInputPerformanceComparison.fig'])
    for fitTypei = 1:length(fitTypes)
        savefig(hModelandData(fitTypei),[saveLocation ...
            '/reducedRankModelVsTargetedDimensions_Ex'...
            num2str(nums(fitTypei)) ...
            '_' fitTypes{fitTypei} '.fig'])
        savefig(gvth(fitTypei) ,[saveLocation ...
            '/reducedRankModelvsGain_Ex'...
            num2str(nums(fitTypei)) ...
            '_' fitTypes{fitTypei} '.fig'])
    end
    
    
end

