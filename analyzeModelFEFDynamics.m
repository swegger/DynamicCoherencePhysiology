function [dK,K, kSteadyState] = analyzeModelFEFDynamics(modelFEF,inputs,varargin)
%%
%
%
%
%%

%% Defaults
plotOpts_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'modelFEF')
addRequired(Parser,'inputs')
addParameter(Parser,'tProbe',[0,150,750])
addParameter(Parser,'kappaSpace',{linspace(-1.5,1.5,48)})
addParameter(Parser,'kappasObserved',NaN)
addParameter(Parser,'plotOpts',plotOpts_default)

parse(Parser,modelFEF,inputs,varargin{:})

modelFEF = Parser.Results.modelFEF;
inputs = Parser.Results.inputs;
tProbe = Parser.Results.tProbe;
kappaSpace = Parser.Results.kappaSpace;
kappasObserved = Parser.Results.kappasObserved;
plotOpts = Parser.Results.plotOpts;

%% Preliminary
N = size(modelFEF.overlaps,1);
M = size(inputs,1);
Sigma = modelFEF.overlaps(:,1:N);
SigmaI = modelFEF.overlaps(:,N+1:end);
sigmas = modelFEF.sigmas;

if length(kappaSpace) == 1
    kappaSpace = repmat(kappaSpace(:),[1,N]);
elseif length(kappaSpace) ~= N
    error(['Space of kappas to explore is ' num2str(length(kappaSpace)) 'D while the model is ' num2str(N) 'D!'])
end

%% Probe the space at selected time points
K = cell(1,N);
[K{:}] = ndgrid(kappaSpace{:});

Klist = nan(numel(K{1}),N);
for i = 1:N
    kappalist(:,i) = K{i}(:);
end

% Preallocate
dkappa = nan(N,size(kappalist,1),length(tProbe));
dK = nan(size(K{1},1),size(K{1},2),length(tProbe));

% Run through each time point and list of kappa values to find the change
% in kappa
for ti = 1:length(tProbe)    
    us = inputs(:,tProbe(ti) == modelFEF.t);
    for listi = 1:size(kappalist,1)
        dkappa(:,listi,ti) = modelDynamics(kappalist(listi,:)',us,Sigma,SigmaI,sigmas);
    end
    dK(:,:,1,ti) = reshape(dkappa(1,:,ti),size(K{1}));
    dK(:,:,2,ti) = reshape(dkappa(2,:,ti),size(K{1}));
end

%% Fixed points

% Trivial fixed point/steady state position
for ti = 1:length(tProbe)
    us = inputs(:,tProbe(ti) == modelFEF.t);
    ssEstimate = kappalist(sum(dkappa(:,:,ti).^2,1) == min(sum(dkappa(:,:,ti).^2,1)),:);
    ssEstimate = ssEstimate(1,:);
    delta = computeDelta(ssEstimate',us,sigmas);
    theta = deltaMean(delta,1e-10,1000);
    SigmaTilde = theta*Sigma;
    SigmaITilde = theta*SigmaI;
    kSteadyState(:,ti) = (eye(N)-SigmaTilde)\SigmaITilde*us;
end

% Nontrivial fixed points
vals = eig(Sigma);
[V,~] = eig(Sigma);
[maxv,maxi] = max(vals);
v = V(:,maxi);
probeLocations = linspace(0.2,1.3,300);
for ti = 1:length(tProbe)
    us = inputs(:,tProbe(ti) == modelFEF.t);
    ssEstimate = kappalist(sum(dkappa(:,:,ti).^2,1) == min(sum(dkappa(:,:,ti).^2,1)),:);
    delta = computeDelta(ssEstimate',us,sigmas);
    for loci = 1:length(probeLocations)
        kappas = kSteadyState(:,ti) + v*probeLocations(loci);
        dkappaTemp(:,loci) = modelDynamics(kappas,us,Sigma,SigmaI,sigmas);
    end
    fp(ti) = probeLocations(sum(dkappaTemp.^2) == min(sum(dkappaTemp.^2)));
end

%% "Null"clines
for ti =1:length(tProbe)
    speed(:,:,ti) = sqrt(dK(:,:,1,ti).^2+dK(:,:,2,ti).^2);
    for k1i = 1:size(K{1},1)
        nullcline2(k1i,:,ti) = [K{1}(k1i,1), K{2}(k1i,speed(k1i,:,ti) == min(speed(k1i,:,ti)))];
    end
    for k2i = 1:size(K{2},2)
        nullcline1(k2i,:,ti) = [K{2}(1,k2i), K{1}(speed(:,k2i,ti) == min(speed(:,k2i,ti)),k2i)];
    end
end

%% Plot results
if plotOpts.On

    hf = figure('Position',[830 127 1174 1062]);
    maxcolumns = 3;
    for ti = 1:length(tProbe)
        if length(tProbe)>maxcolumns
            subplot(ceil(length(tProbe)/maxcolumns),maxcolumns,ti)
        else
            subplot(1,length(tProbe),ti)
        end
        imagesc(kappaSpace{1},kappaSpace{2},log10(speed(:,:,ti)'))
        set(gca,'YDir','normal')
%         h = pcolor(K{1},K{2},log10(speed(:,:,ti)));
%         h.EdgeColor = 'none';
%         h.FaceColor = 'interp';
        caxis([min(log10(speed(:))),max(log10(speed(:)))])
        hold on
        if isfield(plotOpts,'quiverDownSample')
            quiver(K{1}(1:plotOpts.quiverDownSample:end,1:plotOpts.quiverDownSample:end),...
                K{2}(1:plotOpts.quiverDownSample:end,1:plotOpts.quiverDownSample:end),...
                dK(1:plotOpts.quiverDownSample:end,1:plotOpts.quiverDownSample:end,1,ti),...
                dK(1:plotOpts.quiverDownSample:end,1:plotOpts.quiverDownSample:end,2,ti),...
                'k');
        else
            quiver(K{1},K{2},dK(:,:,1,ti),dK(:,:,2,ti),'k');
        end
        plot(kSteadyState(1,ti),kSteadyState(2,ti),'ko','MarkerFaceColor',[1 1 1])
        plot(kSteadyState(1,ti)+v(1)*fp(ti),kSteadyState(2,ti)+v(2)*fp(ti),'ko','MarkerFaceColor',[0 0 0])
        plot(nullcline2(:,1,ti),nullcline2(:,2,ti),'k--')
        plot(nullcline1(:,2,ti),nullcline1(:,1,ti),'k')
        if ~any(isnan(kappasObserved(:)))
            plot(kappasObserved(1,modelFEF.t <= tProbe(ti)),kappasObserved(2,modelFEF.t <= tProbe(ti)),'ro')
        end
        %     plot(kSteadyState(1,ti)+V(1,2),kSteadyState(2,ti)+V(2,2),'ko','MarkerFaceColor',[0 0 0])
        axis equal
        axis([min(K{2}(1,:)) max(K{2}(1,:)) min(K{1}(:,1)) max(K{1}(:,1))])
        title(['t = ' num2str(tProbe(ti))])
    end
    
    if isfield(plotOpts,'saveFile') && isfield(plotOpts,'saveDir')
        
        if ~exist(plotOpts.saveDir,'dir')
            mkdir(plotOpts.saveDir)
        end
        
        savefig(hf,[plotOpts.saveDir '/' plotOpts.saveFile])
        
    end
    
end

%% Functions
%% modelDynamics
function dkappas = modelDynamics(kappas,us,Sigma,SigmaI,sigmas)
    %%
    delta = computeDelta(kappas,us,sigmas);
    theta = deltaMean(delta,1e-10,1000);
    dkappas = -kappas + theta*(Sigma*kappas + SigmaI*us);
    
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
        mu = 1;
    end    
    if any(mu(:) > 1)
        error('Check here')
    end
    
%% tertiaryDelta
function out = tertiaryDelta(deltas,eps,nsteps)
    %%
    zs = linspace(atanh(-1+eps)/deltas,atanh(1-eps)/deltas,nsteps);
    dz = zs(2)-zs(1);
    mu = (1/sqrt(2*pi)) * trapz( dz * exp(-zs.^2/2) .* ...
        (2*(cosh(2*deltas*zs)-2)./cosh(deltas*zs).^4) );