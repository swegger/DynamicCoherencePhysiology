function Beta = MTtoFEFregression(dcp,varargin)
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
trainCondition_default = [true, false, true;
                          true, false, true;
                          true, false, true];
compToBehavioralGain_default.On = false;
theoretical_default.weigthTheory = 'simple';
theoretical_default.expansionDef = 'bestfit';

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'trainCondition',trainCondition_default)
addParameter(Parser,'tWin',[0 900])
addParameter(Parser,'rankN',80)
addParameter(Parser,'ridgeLambda',logspace(-1,12,10))
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'checkMTFit',false)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'zMeanWin',[-Inf,Inf])
addParameter(Parser,'zSTDwin',[-Inf,Inf])
addParameter(Parser,'compToBehavioralGain',compToBehavioralGain_default)
addParameter(Parser,'theoretical',theoretical_default)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
trainCondition = Parser.Results.trainCondition;
directions = Parser.Results.directions;
tWin = Parser.Results.tWin;
rankN = Parser.Results.rankN;
ridgeLambda = Parser.Results.ridgeLambda;
sprefFromFit = Parser.Results.sprefFromFit;
checkMTFit = Parser.Results.checkMTFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
zMeanWin = Parser.Results.zMeanWin;
zSTDwin = Parser.Results.zSTDwin;
compToBehavioralGain = Parser.Results.compToBehavioralGain;
theoretical = Parser.Results.theoretical;
plotOpts = Parser.Results.plotOpts;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Get firing rate data
passCutoff = nan(1000,1);
Rinit = nan(1701,3,3,1000);
cellID = nan(1000,100,3);
indx = 1;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
    
    % Add probe info
    dcp{filei} = addProbeInfo(dcp{filei});
    
    % InitCoh data
    load([sourceDirectory '/' dcp{filei}.datapath(end-8:end-1) 'obj/initCoh' ...
        dcp{filei}.datapath(end-8:end)])
    
    
    if ~isempty(initCoh.R)
        
        passCutoff(indx:indx+length(initCoh.passCutoff)-1) = initCoh.passCutoff;
        
        % Get data for each neuron
        for uniti = 1:length(initCoh.preferredDirectionRelative)
            ind = find(directions == initCoh.preferredDirectionRelative(uniti));
            Rinit(:,:,:,indx) = initCoh.R(:,:,:,uniti,ind);
            
            
            for j = 1:length(initCoh.unitIndex)
                cellID(indx,j,1) = filei;
                cellID(indx,j,2) = initCoh.unitIndex(uniti);
                cellID(indx,j,3) = initCoh.unitIndex(j);
            end
                        
            indx = indx+1;
        end
    end
end


Rinit = Rinit(:,:,:,1:indx-1);
passCutoff = logical(passCutoff(1:indx-1));
cellID = cellID(1:indx-1,:,:);

taccept = initCoh.neuron_t >= tWin(1) & initCoh.neuron_t <= tWin(2);

Rinit = Rinit(taccept,:,:,:);

fef_t = initCoh.neuron_t(taccept);

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff);
cellID = cellID(passCutoff,:,:);


%% Remove outlier rates
m = squeeze(max(Rinit,[],[1,2,3]))*1000;
Rinit = Rinit(:,:,:,m<=150);
cellID = cellID(m<=150,:,:);

%% Get MT data and organize as an input

% Load mtObjs
mtResults = load(objectFile);

% Get mean of each unit

[MT, spref, swidth, dirMT] = getMTdata(mtResults.mt,...
    'speedsMT',speeds,'cohsMT',cohs,'directionsMT',directionsMT,...
    'opponentMT',opponentMT,'sprefFromFit',sprefFromFit,'checkMTFit',checkMTFit,...
    'speedPrefOpts',speedPrefOpts);

% MT = [];
% spref = [];
% for filei = 1:length(mt)
%     disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
%     
%     MTtemp = nan(length(mt{filei}.neuron_t),length(speeds),length(cohs));
%     
%     
%     for si = 1:length(speeds)
%         for ci = 1:length(cohs)
%             [~,condLogical] = trialSort(mt{filei},0,speeds(si),NaN,cohs(ci));
%             MTtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
%         end
%     end
%     
%     if sprefFromFit
%         [mu,sig,~,normR,~] = fitSpeedTuning(mt{filei},'P0',speedPrefOpts.P0,...
%             'ub',speedPrefOpts.ub,'lb',speedPrefOpts.lb,...
%             'c',speedPrefOpts.c,'s',speedPrefOpts.s,'d',speedPrefOpts.d,...
%             'tWin',speedPrefOpts.tWin);
%         spref = [spref mu];
%         
%         if checkMTFit
%             s = linspace(min(speeds)-0.1*min(speeds),max(speeds)*1.1,20);
%             h = figure;
%             semilogx(mt{filei}.speeds(any(mt{filei}.coh(:) == speedPrefOpts.c(:)',2) & any(mt{filei}.directions(:) == speedPrefOpts.d(:)',2)),normR,'o')
%             hold on
%             semilogx(s,speedTuning(mt{filei},s,mu,sig))
%             ax = axis;
%             text(ax(1),0.95*ax(4),['\mu = ' num2str(mu) ', \sig = ' num2str(sig)])
%             input('Press enter to continue ')
%             close(h);
%         end
%     else
%         spref = [spref speeds(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
%             max(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100))))];
%     end
%     MT = cat(4,MT,MTtemp);
% end

mtNeuron_t = mtResults.mt{filei}.neuron_t;


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);

inputs = permute(interpolatedR,[4,1,2,3]);
inputs = inputs(prod(any(isnan(inputs),2),[3,4])==0,:,:,:);

% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,mtNeuron_t);
RtoFit = Rinit(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);
t = mtNeuron_t(iMT);
inputs_z = (inputs-mean(inputs(:,t>zMeanWin(1) & t<=zMeanWin(2),:,:,:),[2,3,4]))./...
    std(inputs(:,t>zSTDwin(1) & t<=zSTDwin(2),:,:),[],[2,3,4]);
RtoFit_z = (RtoFit-mean(RtoFit,[1,2,3]))./std(RtoFit,[],[1,2,3]);

%% Regression models

% Build matrices
Xfit = nan(size(inputs,2)*sum(trainCondition(:)),size(inputs,1));
Yfit = nan(size(RtoFit,1)*sum(trainCondition(:)),size(RtoFit,4));
conditions = nan(size(RtoFit,1)*sum(trainCondition(:)),3);

idx = [1, 1];
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        if trainCondition(si,ci)
            Xfit(idx(1):idx(1)+size(inputs,2)-1,:) = permute(inputs_z(:,:,si,ci),[2,1]);
            Yfit(idx(2):idx(2)+size(RtoFit,1)-1,:) = squeeze(RtoFit_z(:,si,ci,:));
            conditions(idx(2):idx(2)+size(RtoFit,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit,1),1]);
            idx(1) = idx(1) + size(inputs,2);
            idx(2) = idx(2) + size(RtoFit,1);
        end
    end
end

Xtest = nan(size(inputs,2)*sum(~trainCondition(:)),size(inputs,1));
Ytest = nan(size(RtoFit,1)*sum(~trainCondition(:)),size(RtoFit,4));
conditionsTest = nan(size(RtoFit,1)*sum(~trainCondition(:)),3);
idx = [1, 1];
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        if ~trainCondition(si,ci)
            Xtest(idx(1):idx(1)+size(inputs,2)-1,:) = permute(inputs_z(:,:,si,2),[2,1]);
            Ytest(idx(2):idx(2)+size(RtoFit,1)-1,:) = squeeze(RtoFit_z(:,si,2,:));
            conditionsTest(idx(2):idx(2)+size(RtoFit,1)-1,:) = repmat([speedsFEF(si) cohsFEF(ci) trainCondition(si,ci)],[size(RtoFit,1),1]);
            idx(1) = idx(1) + size(inputs,2);
            idx(2) = idx(2) + size(RtoFit,1);
        end
    end
end

% Regress
trainBeta = nan([size(Xfit,2),size(Yfit,2),length(ridgeLambda)]);
for ri = 1:length(ridgeLambda)
    trainBeta(:,:,ri) = inv(Xfit'*Xfit + ridgeLambda(ri)*eye(size(Xfit'*Xfit)))*Xfit'*Yfit;
    trainYfit = Xfit*trainBeta(:,:,ri);
    testYfit = Xtest*trainBeta(:,:,ri);
    
    fitQ(ri,1) = sum( (Yfit(:) - trainYfit(:)).^2 );
    fitQ(ri,2) = sum( (Ytest(:) - testYfit(:)).^2 );
end
Beta = trainBeta(:,:,fitQ(:,2)==min(fitQ(:,2)));

% Predicted responses
YhatTrain = Xfit*Beta;
YhatTest = Xtest*Beta;

% Reduced rank regression
[Uhat,Shat,Vhat] = svd(YhatTrain,0);
BetaRRR = Beta*Vhat(:,1:rankN)*Vhat(:,1:rankN)';
BetaNull = Beta*Vhat(:,end)*Vhat(:,end)';
YrrrTrain = Xfit*BetaRRR;
YnullTrain = Xfit*BetaNull;
YrrrTest = Xtest*BetaRRR;
YnullTest = Xtest*BetaNull;

Rhat = nan(size(RtoFit));
Rnull = nan(size(RtoFit));
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        if trainCondition(si,ci)
            temp = ismember(conditions,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
            Rhat(:,si,ci,:) = YrrrTrain(temp,:);
            Rnull(:,si,ci,:) = YnullTrain(temp,:);
        else
            temp = ismember(conditionsTest,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
            Rhat(:,si,ci,:) = YrrrTest(temp,:);
            Rnull(:,si,ci,:) = YnullTest(temp,:);
        end
    end
end

% Find correlation with held out data and do Fisher transformation
rvals = corrcoef([Ytest YrrrTest]);
RhatCC = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZ = 0.5*log( (1+RhatCC)./(1-RhatCC) );
[~,zvalSort] = sort(RhatZ);

A = Beta*Vhat;
[~,sortInd] = sort(spref(~isnan(interpolatedR(1,1,1,:))));
sp = spref(~isnan(interpolatedR(1,1,1,:)));
sw = swidth(~isnan(interpolatedR(1,1,1,:)));
spref2 = spref(~isnan(interpolatedR(1,1,1,:)));
spref2 = spref2(sortInd);
A = A(sortInd,:);

%% A theoretical treatment, rather than data driven

switch theoretical.weightTheory
    case 'simple'
        % Simple, log2(spref) weighting
        Atheory = [(log2(sp)'-mean(log2(speedsFEF))) ones(size(sp'))];
        
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
        
        Atheory = [uOpt', ones(size(uOpt'))];
end

% Normalize 
Atheory = Atheory./vecnorm(Atheory);

% Construct MT to FEF weight matrix
switch theoretical.expansionDef
    case 'bestfit'
        vOpt = inv(Atheory'*Atheory)*Atheory'*BetaRRR;
        for ai = 1:size(Atheory,2)
            Wopt(:,:,ai) = Atheory(:,ai)*vOpt(ai,:);
        end
        BetaTheory = sum(Wopt,3);
    case 'data'
        for ai = 1:size(Atheory,2)
            for ri = 1 %:rankN
                tempB(:,:,ri) = Atheory(:,ai)*Vhat(:,ri)';
            end
            B(:,:,ai) = sum(tempB,3);
        end
        BetaTheory = sum(B,3);
end
YtheoryTrain = Xfit*BetaTheory;
YtheoryTest = Xtest*BetaTheory;

Rtheory = nan(size(RtoFit));
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        if trainCondition(si,ci)
            temp = ismember(conditions,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
            Rtheory(:,si,ci,:) = YtheoryTrain(temp,:);
        else
            temp = ismember(conditionsTest,[speedsFEF(si),cohsFEF(ci),trainCondition(si,ci)],'rows');
            Rtheory(:,si,ci,:) = YtheoryTest(temp,:);
        end
    end
end
rvals = corrcoef([Ytest YtheoryTest]);
RhatCCTheory = diag(rvals(size(Ytest,2)+1:end,1:size(Ytest,2)));
RhatZTheory = 0.5*log( (1+RhatCCTheory)./(1-RhatCCTheory) );

% Compute output along theoretical channels
for ri = 1:size(Atheory,2)
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            theoreticalOutputs(:,si,ci,ri) = inputs_z(:,:,si,ci)'*Atheory(:,ri);
        end
    end
end

%% Analysis of variance
% Signal dependent noise model: Gaussian of the form exp(-w^2/(2*m*ln(s^2))/sqrt(2*pi*m*ln(s^2))
% rankN = 5;
mhat = sqrt(sum(A(spref2 >= 1,1:rankN).^2./log(spref2(spref2>=1)'.^2),[1,2])./length(spref2(spref2>=1))); % ML estimator of m 
LLsignal_dependent_noise = sum( -log( sqrt(2*pi*log(spref2(spref2>=1)'.^2)) ) -...
    log(mhat) - A(spref2>=1,1:rankN).^2./(2*mhat^2*log(spref2(spref2>=1)'.^2)) ,[1,2]);


% Standard Gaussian noise
sigWeights = sqrt(mean(A(spref2>=1,1:rankN).^2,[1,2]));
LLstandard_noise = sum( -log( sqrt(2*pi) ) -...
    log(sigWeights) - A(spref2>=1,1:rankN).^2./(2*sigWeights^2) ,[1,2]);

%% Plotting
if plotOpts.On

    %% Reduced rank regression analysis
    figure
    for ci = 1:length(cohsFEF)
        for si = 1:length(speedsFEF)
            subplot(length(cohsFEF),length(cohsFEF),si+(ci-1)*length(cohsFEF))
            plot(squeeze(RtoFit_z(:,si,ci,:)),squeeze(Rhat(:,si,ci,:)),'o','Color',speedColors(si,:))
            hold on
            plotUnity;
            xlabel('FEF data (z-score)')
            ylabel('RRR prediction (z-score)')
        end
        
    end
    
    figure
    plot(RhatZ(zvalSort),'o')
    hold on
    plot(RhatZTheory(zvalSort),'o')
    plotHorizontal(2.34/sqrt(size(Ytest,1)-size(Xfit,2)));       % Approximate pval of 0.01
    xlabel('FEF neuron #')
    ylabel('Fisher transformed r-values')
    
    %%
    figure('Name','MT input channels found by RRR','Position',[674 218 1837 1104])
    if rankN < 5
        rankMax = rankN;
    else
        rankMax = 5;
    end
    for ri = 1:rankMax
        for ci = 1:length(cohsFEF)
            subplot(rankMax,length(cohsFEF),ci + (ri-1)*length(cohsFEF))
            for si = 1:length(speedsFEF)
                MToutputChan(:,si,ci,ri) = inputs_z(:,:,si,ci)'*Beta*Vhat(:,ri);
                plot(t,MToutputChan(:,si,ci,ri),'Color',speedColors(si,:))
                hold on
            end
            xlabel('Time from motion onset (ms)')
            ylabel(['Input along dimension#' num2str(ri)])
            tempLims(ci,:) = ylim;
        end
        
        for ci = 1:length(cohsFEF)
            subplot(rankMax,length(cohsFEF),ci + (ri-1)*length(cohsFEF))
            ylim([min(tempLims(:,1)) max(tempLims(:,2))])
        end
    end
    
    %% Weights along 1st reduced rank dimension as a function of preferred speed
    figure
    ranks = 1:4;
    for ri = 1:length(ranks)
        subplot(length(ranks),2,1+(ri-1)*2)
        semilogx(spref2(spref2>=1),A(spref2>=1,ranks(ri)),'o')
        hold on
        semilogx(sp(spref2>=1),Atheory(spref2>=1,1),'o')
        s = logspace(log10(speedPrefOpts.lb(1)),log10(speedPrefOpts.ub(1)),30);
        if LLsignal_dependent_noise > LLstandard_noise
            plot(s,mhat*sqrt(log(s.^2)),'r')
            plot(s,-mhat*sqrt(log(s.^2)),'r')
        end
        
        plot(s,sigWeights*ones(length(s),1),'k')
        plot(s,-sigWeights*ones(length(s),1),'k')
        
        xlim([speedPrefOpts.lb(1), speedPrefOpts.ub(1)])
        ax = axis;
        text(2,0.95*ax(4),['\delta LL = ' num2str(LLsignal_dependent_noise-LLstandard_noise)])
        
        xlabel('Preferred speed (deg/s)')
        ylabel('Weight')
        
        subplot(length(ranks),2,2+(ri-1)*2)
        plot(Vhat(zvalSort,ranks(ri)),'o')
    end
    
    %% Compare output channels to behavioral gain
    if compToBehavioralGain.On
        temp = load(compToBehavioralGain.file,'initGain');
        initGain = temp.initGain;
        
        figure('Name','MT output vs. gain')
        if rankN < 5
            rankMax = rankN;
        else
            rankMax = 5;
        end
        for ri = 1:rankMax
                subplot(1,rankMax,ri)
                for si = 1:length(speedsFEF)                    
                    outputTemp = squeeze(nanmean(MToutputChan(t >= 700 & t <= 800,si,:,ri),1));
                    plot(squeeze(initGain(si,:,3)),outputTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                    hold on
                    outputTemp = squeeze(nanmean(MToutputChan(t >= 100 & t <= 200,si,:,ri),1));
                    plot(squeeze(initGain(si,:,2)),outputTemp,...
                        'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                end
                xlabel('Behavioral gain (unitless)')
                ylabel(['MT input along dimension ' num2str(ri)])
        end
        
        
        figure('Name','Theoretical MT output vs. gain')
        for ri = 1:size(Atheory,2)
                subplot(1,size(Atheory,2),ri)
                for si = 1:length(speedsFEF)                    
                    outputTemp = squeeze(nanmean(theoreticalOutputs(t >= 700 & t <= 800,si,:,ri),1));
                    plot(squeeze(initGain(si,:,3)),outputTemp,...
                        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                    hold on
                    outputTemp = squeeze(nanmean(theoreticalOutputs(t >= 100 & t <= 200,si,:,ri),1));
                    plot(squeeze(initGain(si,:,2)),outputTemp,...
                        'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                end
                xlabel('Behavioral gain (unitless)')
                ylabel(['MT input along dimension ' num2str(ri)])
        end
    end
    
    %% Plot reconstructed firing rates projection on gain dimension
    if compToBehavioralGain.On
        temp = load(compToBehavioralGain.file,'BinitOrth');
        BinitOrth = temp.BinitOrth;
        sz = size(Rhat);
        pInit = reshape((BinitOrth(:,1:4)'*reshape(Rhat,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
        pInitTheory = reshape((BinitOrth(:,1:4)'*reshape(Rtheory,[prod(sz(1:3)),prod(sz(end))])')',[sz(1),sz(2),sz(3),4]);
        
        figure('Name','Predicted response along targeted dimensions','Position',[63 169 1606 1079])
        subplot(3,1,1)
        for speedi = 1:length(speedsFEF)
            for cohi = 1:length(cohsFEF)
                plot(t,pInit(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel('Speed related activity (a.u.)')
        
        subplot(3,1,2)
        for cohi = 1:length(cohsFEF)
            for speedi = 1:length(speedsFEF)
                plot(t,pInit(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel('Coherence related activity (a.u.)')
        
        subplot(3,1,3)
        for cohi = 1:length(cohsFEF)
            for speedi = 1:length(speedsFEF)
                plot(t,pInit(:,speedi,cohi,3),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel('Gain related activity (a.u.)')
        
        
        dimNames = {'Speed','Coherence','Gain','Offset'};
        gvth = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
        for targetedDim = 1:3
            tempData = [];
            subplot(3,1,targetedDim)
            for si = 1:length(speedsFEF)
                initRatesTemp = squeeze(nanmean(pInit(t >= 700 & t <= 800,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                    'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                hold on
                tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
                initRatesTemp = squeeze(nanmean(pInit(t >= 100 & t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
            gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
            axis square
            ax = axis;
            text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
            xlabel('Behavioral gain (unitless)')
            ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
        end
        
        figure('Name','Theory-based response along targeted dimensions','Position',[63 169 1606 1079])
        subplot(3,1,1)
        for speedi = 1:length(speedsFEF)
            for cohi = 1:length(cohsFEF)
                plot(t,pInitTheory(:,speedi,cohi,1),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel('Speed related activity (a.u.)')
        
        subplot(3,1,2)
        for cohi = 1:length(cohsFEF)
            for speedi = 1:length(speedsFEF)
                plot(t,pInitTheory(:,speedi,cohi,2),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel('Coherence related activity (a.u.)')
        
        subplot(3,1,3)
        for cohi = 1:length(cohsFEF)
            for speedi = 1:length(speedsFEF)
                plot(t,pInitTheory(:,speedi,cohi,3),'-','Color',[speedColors(speedi,:) cohsFEF(cohi)/100])
                hold on
            end
        end
        axis tight
        ax(1,:) = axis;
        xlabel('Time from motion onset (ms)')
        ylabel('Gain related activity (a.u.)')
        
        
        dimNames = {'Speed','Coherence','Gain','Offset'};
        gvth2 = figure('Name',['Behavioral gain vs activity on targeted dimension (initCoh)'],'Position',[1956 59 570 1263]);
        for targetedDim = 1:3
            tempData = [];
            subplot(3,1,targetedDim)
            for si = 1:length(speedsFEF)
                initRatesTemp = squeeze(nanmean(pInitTheory(t >= 700 & t <= 800,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,3)),initRatesTemp,...
                    'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                hold on
                tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
                initRatesTemp = squeeze(nanmean(pInitTheory(t >= 100 & t <= 200,si,:,targetedDim),1));
                plot(squeeze(initGain(si,:,2)),initRatesTemp,...
                    'd-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:))
                tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
            end
            gain_projection_cc = corrcoef(tempData(:,1),tempData(:,2));
            axis square
            ax = axis;
            text(0.9*ax(2),0.9*ax(4),['R^2 = ' num2str(gain_projection_cc(1,2).^2)])
            xlabel('Behavioral gain (unitless)')
            ylabel(['Projection on ' dimNames{targetedDim} ' dimension'])
        end
    end
    
    %% Relationship to Functional topography
    if compToBehavioralGain.On
        temp = load(compToBehavioralGain.file,'NumClusters','Y','idx','centInd');
        NumClusters = temp.NumClusters;
        for typei = 1:NumClusters
            typeAngle = wrapTo2Pi(atan2(temp.Y(temp.centInd == typei,2)-mean(temp.Y(:,2)),temp.Y(temp.centInd == typei,1)-mean(temp.Y(:,1))))/2+pi;
            hue = typeAngle / (2 * pi);
            saturation = ones(size(hue));
            value = ones(size(hue));
            hsv = cat(3, hue, saturation, value);
            colorWheel(typei,:) = hsv2rgb(hsv);
        end
        colorWheel(colorWheel > 1) = 1;
        colorWheel(colorWheel < 0) = 0;
        
        
        figure
        ranks = 1:4;
        neuronN = 1:size(Vhat,1);
        for ri = 1:length(ranks)
            subplot(1,length(ranks),ri)
            for typei = 1:NumClusters
                plot(neuronN(temp.idx == typei),Vhat(temp.idx == typei,ri),'o',...
                    'Color',colorWheel(typei,:),'MarkerFaceColor',colorWheel(typei,:))
                hold on
            end
            xlabel('FEF neuron number')
            ylabel(['Loading in dim ' num2str(ri)])
        end
        
        figure
        subplot(1,2,1)
        scatter(temp.Y(:,1),temp.Y(:,2),50,RhatZ,'filled')
        cax(1,:) = caxis;
        xlabel('tSNE_1')
        ylabel('tSNE_2')
        
        subplot(1,2,2)
        scatter(temp.Y(:,1),temp.Y(:,2),50,RhatZTheory,'filled')
        cax(2,:) = caxis;
        xlabel('tSNE_1')
        ylabel('tSNE_2')
        
        subplot(1,2,1)
        caxis([min(cax(:,1)) max(cax(:,2))])
        subplot(1,2,2)
        caxis([min(cax(:,1)) max(cax(:,2))])
        
        figure
        for typei = 1:NumClusters
            subplot(3,NumClusters,typei)
            for si = 1:length(speedsFEF)
                for ci = 1:length(cohsFEF)
                    plot(t,nanmean(RtoFit(:,si,ci,temp.idx==typei),4)*1000,...
                        'Color',[speedColors(si,:) cohsFEF(ci)/100])
                    hold on
                end
            end
            xlabel('Time from motion onset (ms)')
            ylabel('sp/s')
            
            subplot(3,NumClusters,NumClusters+typei)
            imagesc(BetaRRR(sortInd,temp.idx==typei))
            ylabel('MT unit (sorted by speed pref)')
            xlabel('FEF unit in cluster')
            
            subplot(3,NumClusters,2*NumClusters+typei)
            Atemp = Beta(:,temp.idx==typei)*Vhat(temp.idx==typei,:);
            semilogx(spref2,Atemp(sortInd,1),'o')
            ylabel('MT unit weight')
            xlabel('speed preference (deg/s)')
        end
        
        
        figure
        for typei = 1:NumClusters
            subplot(3,NumClusters,typei)
            for si = 1:length(speedsFEF)
                for ci = 1:length(cohsFEF)
                    plot(t,nanmean(Rtheory(:,si,ci,temp.idx==typei),4)*1000,...
                        'Color',[speedColors(si,:) cohsFEF(ci)/100])
                    hold on
                end
            end
            xlabel('Time from motion onset (ms)')
            ylabel('sp/s')
            
            subplot(3,NumClusters,NumClusters+typei)
            imagesc(BetaTheory(sortInd,temp.idx==typei))
            ylabel('MT unit (sorted by speed pref)')
            xlabel('FEF unit in cluster')
            
            subplot(3,NumClusters,2*NumClusters+typei)
            Atemp = BetaTheory(:,temp.idx==typei)*Vhat(temp.idx==typei,:);
            semilogx(spref2,Atemp(sortInd,1:size(Atheory,2)),'o')
            ylabel('MT unit weight')
            xlabel('speed preference (deg/s)')
        end
    end
end