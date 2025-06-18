function compareReducedRankDynamicalSystems(varargin)
%%
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'fitTypes',{'fitSigmas','fitSigmasWithDummy','fitSigmasExtraInput'})
addParameter(Parser,'generationDate',[])

parse(Parser,varargin{:})

fitTypes = Parser.Results.fitTypes;
generationDate = Parser.Results.generationDate;

%% Load results
model_cc.eyeSpeed = nan(length(fitTypes),1000);
model_cc.gain = nan(length(fitTypes),1000);
for fitTypei = 1:length(fitTypes)
    
    sourcePath = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/ar/fitReducedRankModel/' fitTypes{fitTypei} '/'];
    if ~isempty(generationDate)
        sourcePattern = ['*_' generationDate '.mat'];
    else
        sourcePattern = [];
    end
    files = dir([sourcePath, sourcePattern]);
        
    for filei = 1:length(files)
        disp([sourcePath, files(filei).name])
        
        temp = load([sourcePath, files(filei).name],'N','speeds','meanEyeSpeed','eye_t','Rhat','t','initGain');
        for targetedDim = 1:temp.N
            tempData = [];
            for si = 1:length(temp.speeds)
                eyeSpeedTemp = squeeze(nanmean(temp.meanEyeSpeed(temp.eye_t >= 750 & temp.eye_t <= 750,:,:),1));
                initRatesTemp = squeeze(nanmean(temp.Rhat(temp.t >= 700 & temp.t <= 800,si,:,targetedDim),1));
                tempData = [tempData; eyeSpeedTemp(si,:)', temp.initGain(si,:,3)', initRatesTemp(:)];
                
                eyeSpeedTemp = squeeze(nanmean(temp.meanEyeSpeed(temp.eye_t >= 150 & temp.eye_t <= 150,:,:),1));
                initRatesTemp = squeeze(nanmean(temp.Rhat(temp.t >= 100 & temp.t <= 200,si,:,targetedDim),1));
                tempData = [tempData; eyeSpeedTemp(si,:)', temp.initGain(si,:,2)', initRatesTemp(:)];
            end
            rvals = corrcoef(tempData);
            num = regexp(files(filei).name,'_\d+_','match');
            num = str2num(num{1}(num{1}~='_'));
            model_cc.eyeSpeed(fitTypei,num) = rvals(1,3);
            model_cc.gain(fitTypei,num) = rvals(2,3);
        end

    end
end

%% Plot results

figure('Name','Histogram of model output, behavioral data predictive power')
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

figure('Name','Speed vs gain represenation for same MT inputs')
for fitTypei = 1:length(fitTypes)
    subplot(1,length(fitTypes),fitTypei)
    plot(model_cc.eyeSpeed(fitTypei,:),model_cc.gain(fitTypei,:),'o')
    xlabel('Speed prediction performance')
    ylabel('Gain prediction performance')
    axis([-1 1 -1 1])
    axis square
end
    


figure('Name','Across model correlatoins for same MT inputs')
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
    
