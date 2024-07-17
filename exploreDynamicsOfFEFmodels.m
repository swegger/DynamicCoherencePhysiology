function exploreDynamicsOfFEFmodels(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'generateNewInput',false)

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
generateNewInput = Parser.Results.generateNewInput;

%%
switch subject
    %%
    case 'fr'
        
        load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/fitReducedRankModel20240613.mat',...
            'modelFEF','theoreticalInput','objectFile','speeds','cohs','directionsMT','opponentMT',...
            'speedsFEF','cohsFEF','simulateMT','equalizeInputsPriorToStimulusOnset','theoretical',...
            'centerData')
        
        if generateNewInput
            simulateMT.t = modelFEF.t;
            [~, ~, theoreticalInputNew, theoreticalInputNull] = MTpop('objectFile',objectFile,...
                'speeds',speeds,'cohs',cohs,'directionsMT',directionsMT,'opponentMT',opponentMT,...
                'speedsFEF',speedsFEF,'cohsFEF',cohsFEF,'simulateMT',simulateMT,...
                'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
                'theoretical',theoretical,'centerData',centerData);
            theoreticalInputNew(3,:,:,:) = theoreticalInput(3,:,:,:);
            theoreticalInputNull(3,:,:,:) = theoreticalInput(3,:,:,:);
        else
            theoreticalInputNew = theoreticalInput;
        end
        
        % Model
        plotOpts.On = true;
        plotOpts.saveDir = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/fr/reducedRankPhasePlaneAnalysis/' datestr(now,'yyyymmdd')];
        tProbe = ceil(linspace(0,400,9));
%         modelFEF.overlaps(2,5) = 0.1;
        for si = 1:3
            for ci = 1:3
                plotOpts.saveFile = ['phasePlane_speed_' num2str(speedsFEF(si)) '_coh_' num2str(cohsFEF(ci))];
                [~,kappas(:,:,si,ci)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                    't',modelFEF.t,'us',theoreticalInputNew(:,:,si,ci),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
                [dK,K] = analyzeModelFEFDynamics(modelFEF,theoreticalInputNew(:,:,si,ci),'tProbe',tProbe,'kappasObserved',kappas(:,:,si,ci),...
                    'plotOpts',plotOpts);
            end
        end
        
        % Response to null direction
        if generateNewInput
            input = theoreticalInputNull(:,:,3,3);
            [~,kappasNull] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                't',modelFEF.t,'us',input,'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
            [dK,K] = analyzeModelFEFDynamics(modelFEF,input,'tProbe',tProbe,'kappasObserved',kappasNull);
        end
        
        
        % Lesioning inputs
%         input = zeros(size(theoreticalInputNew(:,:,3,3)));
        input = theoreticalInputNew(:,:,1,1);
        input(3,:) = 0;
        [~,kappasZeros] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,'us',input,'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
        [dK,K] = analyzeModelFEFDynamics(modelFEF,input,'tProbe',tProbe,'kappasObserved',kappasZeros);
        
        % Noise in inputs
        for triali = 1:100
            input = theoreticalInputNew(:,:,3,3) + randn(size(theoreticalInputNew)).*sqrt(abs(theoreticalInputNew(:,:,3,1)));
            [~,kappasNoise(:,:,triali)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                't',modelFEF.t,'us',input,'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
        end
        
        

    %%
    case 'ar'
        load('/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240613.mat',...
            'modelFEF','theoreticalInput','objectFile','speeds','cohs','directionsMT','opponentMT',...
            'speedsFEF','cohsFEF','simulateMT','equalizeInputsPriorToStimulusOnset','theoretical',...
            'centerData')
        
        if generateNewInput
            simulateMT.t = modelFEF.t;
            [~, ~, theoreticalInputNew, theoreticalInputNull] = MTpop('objectFile',objectFile,...
                'speeds',speeds,'cohs',cohs,'directionsMT',directionsMT,'opponentMT',opponentMT,...
                'speedsFEF',speedsFEF,'cohsFEF',cohsFEF,'simulateMT',simulateMT,...
                'equalizeInputsPriorToStimulusOnset',equalizeInputsPriorToStimulusOnset,...
                'theoretical',theoretical,'centerData',centerData);
            theoreticalInputNew(3,:,:,:) = theoreticalInput(3,:,:,:);
            theoreticalInputNull(3,:,:,:) = theoreticalInput(3,:,:,:);
        else
            theoreticalInputNew = theoreticalInput;
        end
        
        % Model
        plotOpts.On = true;
        plotOpts.quiverDownSample = 2;
        plotOpts.saveDir = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/ar/reducedRankPhasePlaneAnalysis/' datestr(now,'yyyymmdd')];
        tProbe = ceil(linspace(0,400,9));
%         modelFEF.overlaps(2,5) = 0.1;
        for si = 1:length(speedsFEF)
            for ci = 1:length(cohsFEF)
                plotOpts.saveFile = ['phasePlane_speed_' num2str(speedsFEF(si)) '_coh_' num2str(cohsFEF(ci))];
                [~,kappas(:,:,si,ci)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                    't',modelFEF.t,'us',theoreticalInputNew(:,:,si,ci),'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
                [dK,K] = analyzeModelFEFDynamics(modelFEF,theoreticalInputNew(:,:,si,ci),'tProbe',tProbe,'kappasObserved',kappas(:,:,si,ci),...
                    'plotOpts',plotOpts);
            end
        end
        
        % Response to null direction
        if generateNewInput
            input = theoreticalInputNull(:,:,3,3);
            [~,kappasNull] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                't',modelFEF.t,'us',input,'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
            [dK,K] = analyzeModelFEFDynamics(modelFEF,input,'tProbe',tProbe,'kappasObserved',kappasNull);
        end
        
        % Lesioning inputs
        input = zeros(size(theoreticalInputNew(:,:,3,3)));
%         input = theoreticalInputNew(:,:,1,1);
%         input(3,:) = 0;
        [~,kappasZeros] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,'us',input,'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
        [dK,K] = analyzeModelFEFDynamics(modelFEF,input,'tProbe',tProbe,'kappasObserved',kappasZeros);
        
        % Noise in inputs
        for triali = 1:100
            input = theoreticalInputNew(:,:,3,3) + randn(size(theoreticalInputNew)).*sqrt(abs(theoreticalInputNew(:,:,3,1)));
            [~,kappasNoise(:,:,triali)] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
                't',modelFEF.t,'us',input,'kappas0',modelFEF.R0,'overlaps',modelFEF.overlaps,'sigmas',modelFEF.sigmas);
        end
end