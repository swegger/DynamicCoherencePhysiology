function [m, ste, sensitivity] = compareCohTuningInitDynCoh(initCoh,dynamicCoh,unit,varargin)
%% compareCohTuningInitDynCoh
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'initCoh')
addRequired(Parser,'dynamicCoh')
addRequired(Parser,'unit')
addParameter(Parser,'windows',[50,150; 650,750; 950,1050])
addParameter(Parser,'dirs',0)
addParameter(Parser,'speeds',10)
addParameter(Parser,'cohs',[20,60,100])
addParameter(Parser,'seqs',1:5)

parse(Parser,initCoh,dynamicCoh,unit,varargin{:})

initCoh = Parser.Results.initCoh;
dynamicCoh = Parser.Results.dynamicCoh;
unit = Parser.Results.unit;
windows = Parser.Results.windows;
dirs = Parser.Results.dirs;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
seqs = Parser.Results.seqs;

%% Find counts in each window
for wini = 1:size(windows,1)
    initCounts = cohConditionedCounts(initCoh,'dirs',dirs,'speeds',speeds,'win',windows(wini,:));
    dynCounts = seqConditionedCounts(dynamicCoh,'dirs',0,'win',windows(wini,:));
    
    % Limit trial number to max w/out NaN
    maxtrials = min(min(min(sum(~isnan(initCounts),1))));
    initCounts = initCounts(1:maxtrials,:,:,:);
    maxtrials = min(min(sum(~isnan(dynCounts),1)));
    dynCounts = dynCounts(1:maxtrials,:,:);
    
    % Find mean and standard error
    m.init(:,wini) = mean(initCounts(:,:,:,initCoh.unitIndex == unit),1);
    ste.init(:,wini) = std(initCounts(:,:,:,initCoh.unitIndex == unit),[],1)/sqrt(size(initCounts,1));
    m.dyn(:,wini) = mean(dynCounts(:,:,dynamicCoh.unitIndex == unit),1);
    ste.dyn(:,wini) = std(dynCounts(:,:,dynamicCoh.unitIndex == unit),[],1)/sqrt(size(dynCounts,1));
    for seqi = 1:length(seqs)
        cohmap(seqi,wini) = dynamicCoh.coh(round(mean(windows(wini,:)))+101,seqi);
    end
    
    % Fit LNP model
    ind = 0;
    for triali = 1:size(initCounts,1)
        for cohi = 1:length(cohs)
            for speedi = 1:length(speeds)
                ind = ind+1;
                X(ind,:) = cohs(cohi);
                c(ind) = initCounts(triali,speedi,cohi,initCoh.unitIndex == unit);
            end
        end
    end
    [theta, logPosterior] = FitLinearExpPoissonObsMAP(X,c',[0 0]);
    [thetaC, logPosteriorC] = FitConstantPoissonObsMAP(X,c');
    
    sensitivity.init.theta(:,wini) = theta;
    sensitivity.init.deltaBIC(:,wini) = 2*(logPosterior - 2*log(size(initCounts,1))/2) - ...
        2*(logPosteriorC - 2*log(size(initCounts,1)));
    
    clear X c
    
    ind = 0;
    for triali = 1:size(dynCounts,1)
        for seqi = 1:length(seqs)
                ind = ind+1;
                X(ind,:) = dynamicCoh.coh(round(mean(windows(wini,:)))+101,seqi);
                c(ind) = dynCounts(triali,seqi,dynamicCoh.unitIndex == unit);
        end
    end
    [theta, logPosterior] = FitLinearExpPoissonObsMAP(X,c',[0 0]);
    [~, logPosteriorC] = FitConstantPoissonObsMAP(X,c');
    
    sensitivity.dyn.theta(:,wini) = theta;
    sensitivity.dyn.deltaBIC(:,wini) = 2*(logPosterior - 2*log(size(dynCounts,1))/2) - ...
        2*(logPosteriorC - 2*log(size(dynCounts,1)));
    
    clear X c
    
end