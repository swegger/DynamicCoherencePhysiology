function lesionRRdynamicsModel(modelFEF,inputs,varargin)
%% lesionRRdynamicsModel
%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'modelFEF')
addRequired(Parser,'inputs')
addParameter(Parser,'speedsFEF',[5 10 20])
addParameter(Parser,'cohsFEF',[20 60 100])

parse(Parser,modelFEF,inputs,varargin{:})

modelFEF = Parser.Results.modelFEF;
inputs = Parser.Results.inputs;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;

%% Preliminary

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speedsFEF),'sampleN',length(speedsFEF));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% No lesion
overlapsTemp = modelFEF.overlaps;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Unlesioned')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%%
figure('Name','Unlesioned fixed points with inputs')
ind = 0;
% [Sgrid,Ggrid] = meshgrid(linspace(0,1.2,10),linspace(-1.2,0,10));
% sgrid = Sgrid(:);
% ggrid = Ggrid(:);
pertSize = 0.2;
tFixed = 650;
for ci = 1:size(inputs,4)
%         ind = ind+1;
%         subplot(size(inputs,4),1,ind)
        for si = 1:size(inputs,3)
            I = eye(size(sigmaTildes,1));
            kI = inputs(:,modelFEF.t == tFixed,si,ci);
            Sigma = sigmaTildes(:,1:2,modelFEF.t==tFixed,si,ci);
            SigmaI = sigmaTildes(:,2+1:end,modelFEF.t==tFixed,si,ci);
            theta = SigmaI(1,1)./modelFEF.overlaps(1,3);
            [A,B] = eig(Sigma/theta);
            kappasSteadyState(:,si,ci) = (I-Sigma)\SigmaI*kI;
            
            [Sgrid,Ggrid] = meshgrid(linspace(kappasSteadyState(1,si,ci)-pertSize,kappasSteadyState(1,si,ci)+pertSize,10),...
                linspace(kappasSteadyState(2,si,ci)-pertSize,kappasSteadyState(2,si,ci)+pertSize,10));
            sgrid = Sgrid(:);
            ggrid = Ggrid(:);
            
            dk(:,:,si,ci) = -[sgrid, ggrid]' + Sigma*[sgrid, ggrid]' + SigmaI*kI;
            quiver(Sgrid,Ggrid,reshape(dk(1,:,si,ci),size(Sgrid)),reshape(dk(2,:,si,ci),size(Sgrid)),...
                'Color',speedColors(si,:))
            hold on
            plot(kappasSteadyState(1,si,ci),kappasSteadyState(2,si,ci),'ko')
            plot(kappas(1,modelFEF.t==tFixed,si,ci),kappas(2,modelFEF.t==tFixed,si,ci),'o','Color',speedColors(si,:))
            plot([kappasSteadyState(1,si,ci),kappas(1,modelFEF.t==tFixed,si,ci)],...
                [kappasSteadyState(2,si,ci),kappas(2,modelFEF.t==tFixed,si,ci)],'k-')
        end
end

%% Lesion input 1
overlapsTemp = modelFEF.overlaps;
inputTemp = inputs;
inputTemp(1,:,:,:) = 0;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion input 1')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Lesion input 2
overlapsTemp = modelFEF.overlaps;
inputTemp = inputs;
inputTemp(2,:,:,:) = 0;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion input 2')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end


%% Lesion input 3
overlapsTemp = modelFEF.overlaps;
inputTemp = inputs;
inputTemp(3,:,:,:) = 0;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion input 3')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end


%% Lesion recurrence 1,1
overlapsTemp = modelFEF.overlaps;
overlapsTemp(1,1) = 0;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion overlaps(1,1)')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Lesion recurrence 1,2
overlapsTemp = modelFEF.overlaps;
overlapsTemp(1,2) = 0;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion overlaps(1,2)')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Lesion recurrence 2,1
overlapsTemp = modelFEF.overlaps;
overlapsTemp(2,1) = 0;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion overlaps(2,1)')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Lesion recurrence 2,2
overlapsTemp = modelFEF.overlaps;
overlapsTemp(2,2) = 0;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion overlaps(2,2)')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Lesion nonlinearity
overlapsTemp = modelFEF.overlaps;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas,...
            'forceLinear',true);
            
    end
end

figure('Name','Lesion nonlinearity')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end


%% Lesion recurrence 1,2 and 2,1
overlapsTemp = modelFEF.overlaps;
overlapsTemp(1,2) = 0;
overlapsTemp(2,1) = 0;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion overlaps(1,2) and overlaps(2,1)')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Lesion recurrence 1,2 and 2,2
overlapsTemp = modelFEF.overlaps;
overlapsTemp(1,2) = 0;
overlapsTemp(2,2) = 0;
inputTemp = inputs;

% overlapsTemp(2,N+1) = 1.7;
for si = 1:length(speedsFEF)
    for ci = 1:length(cohsFEF)
        [~,kappas(:,:,si,ci), sigmaTildes(:,:,:,si,ci), deltas(:,si,ci),~] = simulateLatentDynamics('tau',modelFEF.tau/modelFEF.dt,...
            't',modelFEF.t,...
            'us',inputTemp(:,:,si,ci),...
            'kappas0',modelFEF.R0,...
            'overlaps',overlapsTemp,...
            'sigmas',modelFEF.sigmas);
            
    end
end

figure('Name','Lesion overlaps(1,2) and overlaps(2,2)')
for ki = 1:size(kappas,1)
    for ci = 1:size(inputs,4)
        subplot(size(kappas,1),size(inputs,4),ci+(ki-1)*size(inputs,4))
        for si = 1:size(inputs,3)
            plot(modelFEF.t,kappas(ki,:,si,ci),...
                'Color',speedColors(si,:))
            hold on
        end
        xlabel('Time from motion onset (ms)')
        ylabel(['Acitivity on ' modelFEF.dimNames{ki}])
        ylim([-1.1 1.1])
        xlim([modelFEF.t(1) modelFEF.t(end)])
    end
end

%% Functions

%% computeDelta
function delta = computeDelta(kappas,us,sigmas)
    %%
    a = sum(kappas.^2 .* sigmas(1:size(kappas,1)).^2) + ...
        sum(us.^2 .* sigmas(size(kappas,1)+1:size(kappas,1)+size(us,1)).^2);
    delta = sqrt(a);
    
%% deltaMean
function mu = deltaMean(deltas,eps,nsteps)
    %%
    if deltas>0
        zs = linspace(atanh(-1+eps)/deltas,atanh(1-eps)/deltas,nsteps);
        dz = zs(2)-zs(1);
        mu = (1/sqrt(2*pi)) * trapz( dz * exp(-zs.^2/2) .* ...
            (1 - tanh(deltas*zs).^2) );
    else
        mu = 0;
    end    
    if any(mu(:) > 1)
        error('Check here')
    end