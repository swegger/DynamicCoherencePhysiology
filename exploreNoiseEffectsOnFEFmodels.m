function exploreNoiseEffectsOnFEFmodels(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'lambdaS',31*0.3)
addParameter(Parser,'lambdaT',100)
addParameter(Parser,'mtSig',logspace(-3,-0.5,4))
addParameter(Parser,'otherInputSig',logspace(-2,0,4))
addParameter(Parser,'trialN',200)
addParameter(Parser,'noiseTimeConstant',20/1000)
addParameter(Parser,'otherInputNoiseModel','lowpassNoise')

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
lambdaS = Parser.Results.lambdaS;
mtSig = Parser.Results.mtSig;
otherInputSig = Parser.Results.otherInputSig;
trialN = Parser.Results.trialN;
noiseTimeConstant = Parser.Results.noiseTimeConstant;
otherInputNoiseModel = Parser.Results.otherInputNoiseModel;

%% Run simulations
switch subject
    %%
    case 'fr'
        load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/fitReducedRankModel20240613.mat',...
            'modelFEF','theoreticalInput','objectFile','speeds','cohs','directionsMT','opponentMT',...
            'speedsFEF','cohsFEF','simulateMT','equalizeInputsPriorToStimulusOnset','theoretical',...
            'centerData','meanEyeSpeed','initGain','eye_t')
        
        simulateMT.t = modelFEF.t;
        [~, inputs_n, ~,~,spref] = MTpop('objectFile',objectFile,...
            'speeds',speeds,'cohs',cohs,'directionsMT',directionsMT,'opponentMT',opponentMT,...
            'speedsFEF',speedsFEF,'cohsFEF',cohsFEF,'simulateMT',simulateMT,...
            'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
            'theoretical',theoretical,'centerData',centerData);
        
        % Noiseless inputs
        Atheory = [(log2(spref)' - log2(mean(speedsFEF))) ones(size(spref'))];
        Atheory = Atheory./vecnorm(Atheory);
        for ci = 1:length(cohsFEF)
            for si = 1:length(speedsFEF)
                for ri = 1:size(Atheory,2)
                    input(ri,:,si,ci) = (inputs_n(:,:,si,ci))'*...
                        Atheory(:,ri);
                end
%                 input(3,:,si,ci) = theoreticalInput(3,:,si,ci);
                input(3,:,si,ci) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2))*modelFEF.extraInput(1) ));
            end
        end
        
        inputNormalized = input./max(abs(input),[],[2,3,4]);
        for si = 1:3
            for ci = 1:3
                [~,kappas(:,:,si,ci)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                    't',modelFEF.t,'us',inputNormalized(:,:,si,ci),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
            end
        end
        
        
        % Simulate the model with noise as specified
        kappasNoise = nan(size(modelFEF.overlaps,1),size(simulateMT.t,2),trialN,length(speedsFEF),length(cohsFEF));
        for mtSigi = 1:length(mtSig)
            for inSigi = 1:length(otherInputSig)
                % Construct covariance matrix
                %         Sigma = covarainceFunctionFull(lambdaS,lambdaT,simulateMT.t,spref);
                Sigma = mtSig(mtSigi)^2*covarainceFunctionSpeedPreference(lambdaS,spref);
                
               % Simulate noise
                for ci = 3%1:length(cohsFEF)
                    for si = 1:length(speedsFEF)
                        noise = mvnrnd(zeros(size(Sigma,1),1),Sigma,trialN*length(simulateMT.t));
                        N = reshape(noise,[length(simulateMT.t) trialN length(spref)]);
                        
                        % Noise in inputs
                        for triali = 1:trialN
                            Nlp = lowpass(squeeze(N(:,triali,:)),1/noiseTimeConstant,1000);
                            
                            for ri = 1:size(Atheory,2)
                                inputNoise(ri,:,si,ci) = (inputs_n(:,:,si,ci)+Nlp')'*...
                                    Atheory(:,ri);
                            end
                            
                            switch otherInputNoiseModel
                                case 'lowpassNoise'
                                    Nlp_otherInput = lowpass(otherInputSig(inSigi)*randn(size(N,1),1),1/noiseTimeConstant,1000);
                                    inputNoise(3,:,si,ci) = theoreticalInput(3,:,si,ci) + Nlp_otherInput';
                                case 'rampSlope'
                                    inputNoise(3,:,si,ci) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2))*modelFEF.extraInput(1)*Nlp_otherInput ));
                                    Nlp_otherInput = gamrnd(1/otherInputSig(inSigi),1)/((1/otherInputSig(inSigi)-1));
                                case 'rampMidpoint'
                                    inputNoise(3,:,si,ci) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2)-Nlp_otherInput)*modelFEF.extraInput(1) ));
                                    Nlp_otherInput = randn*otherInputSig(inSigi);
                                otherwise
                                    error(['otherInputNoiseModel ' otherInputNoiseModel ' not recognized!'])
                            end
                            inputNoise(:,:,si,ci) = inputNoise(:,:,si,ci)./max(abs(input),[],[2,3,4]);
                            
                            [~,kappasNoise(:,:,triali,si,ci,mtSigi,inSigi)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                                't',modelFEF.t,'us',inputNoise(:,:,si,ci),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
                            
                        end
                    end
                end
                
                speedNoise(mtSigi,inSigi,1) = var(kappasNoise(1,simulateMT.t==150,:,:,3,mtSigi,inSigi),[],'all');
                gainNoise(mtSigi,inSigi,1) = var(kappasNoise(2,simulateMT.t == 150,:,:,3,mtSigi,inSigi),[],'all');
                speedNoise(mtSigi,inSigi,2) = var(kappasNoise(1,simulateMT.t==750,:,:,3,mtSigi,inSigi),[],'all');
                gainNoise(mtSigi,inSigi,2) = var(kappasNoise(2,simulateMT.t == 750,:,:,3,mtSigi,inSigi),[],'all');
%                 
%                 figure
%                 subplot(1,2,1)
%                 plot(simulateMT.t,squeeze(kappas(1,:,:,3)))
%                 hold on
%                 plot(simulateMT.t,squeeze(kappasNoise(1,:,:,2,3,mtSigi,inSigi)),'k')
%                 
%                 subplot(1,2,2)
%                 plot(simulateMT.t,squeeze(kappas(2,:,:,3)))
%                 hold on
%                 plot(simulateMT.t,squeeze(kappasNoise(2,:,:,2,3,mtSigi,inSigi)),'k')
                
            end
        end
        
        % Build regression model to relate kappa to speed and gain
        t = simulateMT.t;
        tempData = [];
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
            initRatesTemp = squeeze(nanmean(kappas(1,t >= 700 & t <= 800,si,:),2));
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
            initRatesTemp = squeeze(nanmean(kappas(1,t >= 100 & t <= 200,si,:),2));
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        speedBeta = regress(tempData(:,1),[tempData(:,2),ones(size(tempData,1),1)]);
        
        tempData = [];
        for si = 1:length(speedsFEF)
            initRatesTemp = squeeze(nanmean(kappas(2,t >= 700 & t <= 800,si,:),2));
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
            initRatesTemp = squeeze(nanmean(kappas(2,t >= 100 & t <= 200,si,:),2));
            tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
        end
        gainBeta = regress(tempData(:,1),[tempData(:,2),ones(size(tempData,1),1)]);
        
        
    case 'ar'
        load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240613.mat',...
            'modelFEF','theoreticalInput','objectFile','speeds','cohs','directionsMT','opponentMT',...
            'speedsFEF','cohsFEF','simulateMT','equalizeInputsPriorToStimulusOnset','theoretical',...
            'centerData','meanEyeSpeed','initGain','eye_t')
        
        simulateMT.t = modelFEF.t;
        [~, inputs_n, ~,~,spref] = MTpop('objectFile',objectFile,...
            'speeds',speeds,'cohs',cohs,'directionsMT',directionsMT,'opponentMT',opponentMT,...
            'speedsFEF',speedsFEF,'cohsFEF',cohsFEF,'simulateMT',simulateMT,...
            'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
            'theoretical',theoretical,'centerData',centerData);
        
        % Noiseless inputs
        Atheory = [(log2(spref)' - log2(mean(speedsFEF))) ones(size(spref'))];
        Atheory = Atheory./vecnorm(Atheory);
        for ci = 1:length(cohsFEF)
            for si = 1:length(speedsFEF)
                for ri = 1:size(Atheory,2)
                    input(ri,:,si,ci) = (inputs_n(:,:,si,ci))'*...
                        Atheory(:,ri);
                end
%                 input(3,:,si,ci) = theoreticalInput(3,:,si,ci);
                input(3,:,si,ci) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2))*modelFEF.extraInput(1) ));
            end
        end
        
        inputNormalized = input./max(abs(input),[],[2,3,4]);
        for si = 1:3
            for ci = 1:3
                [~,kappas(:,:,si,ci)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                    't',modelFEF.t,'us',inputNormalized(:,:,si,ci),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
            end
        end
        
        
        % Simulate the model with noise as specified
        kappasNoise = nan(size(modelFEF.overlaps,1),size(simulateMT.t,2),trialN,length(speedsFEF),length(cohsFEF));
        for mtSigi = 1:length(mtSig)
            for inSigi = 1:length(otherInputSig)
                % Construct covariance matrix
                %         Sigma = covarainceFunctionFull(lambdaS,lambdaT,simulateMT.t,spref);
                Sigma = mtSig(mtSigi)^2*covarainceFunctionSpeedPreference(lambdaS,spref);
                
               % Simulate noise
                for ci = 3%1:length(cohsFEF)
                    for si = 1:length(speedsFEF)
                        noise = mvnrnd(zeros(size(Sigma,1),1),Sigma,trialN*length(simulateMT.t));
                        N = reshape(noise,[length(simulateMT.t) trialN length(spref)]);
                        
                        % Noise in inputs
                        for triali = 1:trialN
                            Nlp = lowpass(squeeze(N(:,triali,:)),1/noiseTimeConstant,1000);
                            
                            for ri = 1:size(Atheory,2)
                                inputNoise(ri,:,si,ci,triali) = (inputs_n(:,:,si,ci)+Nlp')'*...
                                    Atheory(:,ri);
                            end
                            
                            switch otherInputNoiseModel
                                case 'lowpassNoise'
                                    Nlp_otherInput = lowpass(otherInputSig(inSigi)*randn(size(N,1),1),1/noiseTimeConstant,1000);
                                    inputNoise(3,:,si,ci,triali) = theoreticalInput(3,:,si,ci) + Nlp_otherInput';
                                case 'rampSlope'
                                    Nlp_otherInput = gamrnd(1/otherInputSig(inSigi),1)/((1/otherInputSig(inSigi)-1));
                                    inputNoise(3,:,si,ci,triali) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2))*modelFEF.extraInput(1)*Nlp_otherInput ));
                                case 'rampMidpoint'
                                    Nlp_otherInput = randn*otherInputSig(inSigi);
                                    inputNoise(3,:,si,ci,triali) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2)-Nlp_otherInput)*modelFEF.extraInput(1) ));
                                case 'rampSlope&Midpoint'
                                    Nlp_otherInput = gamrnd(1/otherInputSig(inSigi),1)/((1/otherInputSig(inSigi)-1));
                                    inputNoise(3,:,si,ci,triali) = 1./(1 + exp( -(simulateMT.t-modelFEF.extraInput(2)-randn*20)*modelFEF.extraInput(1)*Nlp_otherInput ));                                    
                                otherwise
                                    error(['otherInputNoiseModel ' otherInputNoiseModel ' not recognized!'])
                            end
                            inputNoise(:,:,si,ci,triali) = inputNoise(:,:,si,ci,triali)./max(abs(input),[],[2,3,4]);
                            
                            [~,kappasNoise(:,:,triali,si,ci,mtSigi,inSigi)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                                't',modelFEF.t,'us',inputNoise(:,:,si,ci,triali),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
                            
                        end
                    end
                end
                
%                 speedNoise(mtSigi,inSigi,1) = var(kappasNoise(1,simulateMT.t==150,:,:,3,mtSigi,inSigi),[],'all');
                kappasTemp = squeeze(kappasNoise(1,simulateMT.t==150,:,:,3,mtSigi,inSigi));
                v = robustSpeedVarianceEstimate(kappasTemp,20);
                speedNoise(mtSigi,inSigi,1) = sum(v);
                
%                 gainNoise(mtSigi,inSigi,1) = var(kappasNoise(2,simulateMT.t == 150,:,:,3,mtSigi,inSigi),[],'all');
                kappasTemp = squeeze(kappasNoise(2,simulateMT.t==150,:,:,3,mtSigi,inSigi));
                v = robustGainVarianceEstimate(kappasTemp,20);
                gainNoise(mtSigi,inSigi,1) = v;
                
%                 speedNoise(mtSigi,inSigi,2) = var(kappasNoise(1,simulateMT.t==750,:,:,3,mtSigi,inSigi),[],'all');
                kappasTemp = squeeze(kappasNoise(1,simulateMT.t==750,:,:,3,mtSigi,inSigi));
                v = robustSpeedVarianceEstimate(kappasTemp,20);
                speedNoise(mtSigi,inSigi,2) = sum(v);
                
%                 gainNoise(mtSigi,inSigi,2) = var(kappasNoise(2,simulateMT.t == 750,:,:,3,mtSigi,inSigi),[],'all');
                kappasTemp = squeeze(kappasNoise(2,simulateMT.t==750,:,:,3,mtSigi,inSigi));
                v = robustGainVarianceEstimate(kappasTemp,20);
                gainNoise(mtSigi,inSigi,2) = v;
%                 
%                 figure
%                 subplot(1,2,1)
%                 plot(simulateMT.t,squeeze(kappas(1,:,:,3)))
%                 hold on
%                 plot(simulateMT.t,squeeze(kappasNoise(1,:,:,2,3,mtSigi,inSigi)),'k')
%                 
%                 subplot(1,2,2)
%                 plot(simulateMT.t,squeeze(kappas(2,:,:,3)))
%                 hold on
%                 plot(simulateMT.t,squeeze(kappasNoise(2,:,:,2,3,mtSigi,inSigi)),'k')
                
            end
        end
        
        % Build regression model to relate kappa to speed and gain
        t = simulateMT.t;
        tempData = [];
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 750 & eye_t <= 750,:,:),1));
            initRatesTemp = squeeze(nanmean(kappas(1,t >= 700 & t <= 800,si,:),2));
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        for si = 1:length(speeds)
            eyeSpeedTemp = squeeze(nanmean(meanEyeSpeed(eye_t >= 150 & eye_t <= 150,:,:),1));
            initRatesTemp = squeeze(nanmean(kappas(1,t >= 100 & t <= 200,si,:),2));
            tempData = [tempData; eyeSpeedTemp(si,:)',initRatesTemp(:)];
        end
        speedBeta = regress(tempData(:,1),[tempData(:,2),ones(size(tempData,1),1)]);
        
        tempData = [];
        for si = 1:length(speedsFEF)
            initRatesTemp = squeeze(nanmean(kappas(2,t >= 700 & t <= 800,si,:),2));
            tempData = [tempData; initGain(si,:,3)',initRatesTemp(:)];
            initRatesTemp = squeeze(nanmean(kappas(2,t >= 100 & t <= 200,si,:),2));
            tempData = [tempData; initGain(si,:,2)',initRatesTemp(:)];
        end
        gainBeta = regress(tempData(:,1),[tempData(:,2),ones(size(tempData,1),1)]);
end

%% Ploting
figure
[INSIG,MTSIG] = meshgrid(otherInputSig,mtSig);
subplot(2,2,1)
contour(MTSIG,INSIG,sqrt(speedNoise(:,:,1)*speedBeta(1)^2),logspace(log10(min(sqrt(speedNoise(:)*speedBeta(1)^2))),log10(max(sqrt(speedNoise(:)*speedBeta(1)^2))),10))
hold on
contour(MTSIG,INSIG,sqrt(speedNoise(:,:,1)*speedBeta(1)^2),sqrt([(0.043^2*sum(speedsFEF.^2)), (0.094^2*sum(speedsFEF.^2))]),'r','LineWidth',2)
xlabel('MT input noise')
ylabel('Extra input noise')
title('Speed noise, early')
set(gca,'XScale','log','YScale','log')
axis square
colorbar

subplot(2,2,2)
contour(MTSIG,INSIG,sqrt(gainNoise(:,:,1)*gainBeta(1)^2),logspace(log10(min(sqrt(gainNoise(:)*gainBeta(1)^2))),log10(max(sqrt(gainNoise(:)*gainBeta(1)^2))),10))
hold on
contour(MTSIG,INSIG,sqrt(gainNoise(:,:,1)*gainBeta(1)^2),[0.043 0.07],'r','LineWidth',2)
xlabel('MT input noise')
ylabel('Extra input noise')
title('Gain noise, early')
set(gca,'XScale','log','YScale','log')
axis square
colorbar

subplot(2,2,3)
contour(MTSIG,INSIG,sqrt(speedNoise(:,:,2)*speedBeta(1)^2),logspace(log10(min(sqrt(speedNoise(:)*speedBeta(1)^2))),log10(max(sqrt(speedNoise(:)*speedBeta(1)^2))),10))
hold on
contour(MTSIG,INSIG,sqrt(speedNoise(:,:,2)*speedBeta(1)^2),sqrt([(0.043^2*sum(speedsFEF.^2)), (0.094^2*sum(speedsFEF.^2))]),'r','LineWidth',2)
xlabel('MT input noise')
ylabel('Extra input noise')
title('Speed noise, late')
set(gca,'XScale','log','YScale','log')
axis square
colorbar

subplot(2,2,4)
contour(MTSIG,INSIG,sqrt(gainNoise(:,:,2)*gainBeta(1)^2),logspace(log10(min(sqrt(gainNoise(:)*gainBeta(1)^2))),log10(max(sqrt(gainNoise(:)*gainBeta(1)^2))),10))
hold on
contour(MTSIG,INSIG,sqrt(gainNoise(:,:,2)*gainBeta(1)^2),[0.043 0.07],'r','LineWidth',2)
xlabel('MT input noise')
ylabel('Extra input noise')
title('Gain noise, late')
set(gca,'XScale','log','YScale','log')
axis square
colorbar

%% functions

% covarianceFunction
function Sigma = covarainceFunctionFull(lambdaS,lambdaT,t,prefSpeeds)
    t = t(:);
    prefSpeeds = prefSpeeds(:)';
    allSpeeds = repmat(prefSpeeds,[1,length(t)]);
    allTimes = repmat(t,[1,length(prefSpeeds)]);
    allTimes = reshape(allTimes',[1,length(allSpeeds)]);
    
    Sigma = exp(-(allSpeeds-allSpeeds').^2/2/lambdaS^2) .* exp(-abs(allTimes-allTimes')/lambdaT);
    
    
function Sigma = covarainceFunctionSpeedPreference(lambdaS,prefSpeeds)
    prefSpeeds = prefSpeeds(:)';
    
    Sigma = exp(-(prefSpeeds-prefSpeeds').^2/2/lambdaS^2);
    
function v = robustSpeedVarianceEstimate(kappasTemp,countMax)
    outN = Inf;
    count = 0;
    while outN > 10 && count < countMax
        count = count+1;
        m = nan(1,size(kappasTemp,2));
        v = nan(1,size(kappasTemp,2));
        kappaOut = false(size(kappasTemp));
        for si = 1:size(kappasTemp,2)
            m(si) = nanmean(kappasTemp(:,si));
            v(si) = nanvar(kappasTemp(:,si));
            kappaOut(:,si) = kappasTemp(:,si) <= (m(si)+2*sqrt(v(si))) & kappasTemp(:,si) >= (m(si)-2*sqrt(v(si)));
        end
        kappasTemp(~kappaOut) = NaN;
        outN = sum(~kappaOut(:));
    end
    
function v = robustGainVarianceEstimate(kappasTemp,countMax)
    outN = Inf;
    count = 0;
    while outN > 10 && count <= countMax
        count = count+1;
        m = mean(kappasTemp,'all');
        v = var(kappasTemp,[],'all');
        kappasTemp = kappasTemp(kappasTemp <= (m+2*sqrt(v)) & kappasTemp >= (m-2*sqrt(v)));
        outN = sum(kappasTemp <= (m+2*sqrt(v)) & kappasTemp >= (m-2*sqrt(v)));
    end
    