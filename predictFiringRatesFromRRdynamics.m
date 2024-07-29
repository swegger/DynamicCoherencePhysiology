function predictFiringRatesFromRRdynamics(subject,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'subject')
addParameter(Parser,'neuronTypingFile',NaN)
addParameter(Parser,'fitReducedRankDynamicsFile',NaN)
addParameter(Parser,'MTtoFEFregressionFile',NaN)
addParameter(Parser,'speedCluster',9)
addParameter(Parser,'speeds',[5,10,20])
addParameter(Parser,'cohs',[20 60 100])
addParameter(Parser,'tol',1)
addParameter(Parser,'ridgeLambda',0)

parse(Parser,subject,varargin{:})

subject = Parser.Results.subject;
neuronTypingFile = Parser.Results.neuronTypingFile;
fitReducedRankDynamicsFile = Parser.Results.fitReducedRankDynamicsFile;
MTtoFEFregressionFile = Parser.Results.MTtoFEFregressionFile;
speedCluster = Parser.Results.speedCluster;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
tol = Parser.Results.tol;
ridgeLambda = Parser.Results.ridgeLambda;

%% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
figure;
colors = colormap('lines');
close(gcf)
initColors = 1-[20 20 20; 60 60 60; 100 100 100]/100;

%% Check for specific files as input, otherwise set default files
if any(isnan(neuronTypingFile))
    switch subject
        case 'ar'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';
        case 'fr'
            neuronTypingFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/arfr/neuronTypingAnalysis/neuronTyping20240530.mat';            
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(fitReducedRankDynamicsFile))
    switch subject
        case 'ar'
%             fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240613.mat';
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/ar/ReducedRankModel/fitReducedRankModel20240729.mat';
        case 'fr'
            fitReducedRankDynamicsFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFdynamics/fr/ReducedRankModel/fitReducedRankModel20240613.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

if any(isnan(MTtoFEFregressionFile))
    switch subject
        case 'ar'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/ar/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        case 'fr'
            MTtoFEFregressionFile = '/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/fr/MTtoFEFregressionResults/MTtoFEFregressionResults20240603.mat';
        otherwise
            error(['Subejct ' subject ' not recognized!'])
    end
end

%% Load the data
neuronTyping = load(neuronTypingFile,'subjects','Y','cellID','gainRegression','idx','Rinit','NumClusters','initCoh');
rrModel = load(fitReducedRankDynamicsFile,'modelFEF','theoreticalInput','kappas');
MTtoFEF = load(MTtoFEFregressionFile,'vOpt');

%% Select data from neuronTyping related to this subject
subjecti = find(strcmp(neuronTyping.subjects,subject));
Y = neuronTyping.Y(neuronTyping.cellID(:,1,4) == subjecti,:);
idx = neuronTyping.idx(neuronTyping.cellID(:,1,4) == subjecti);
R = neuronTyping.Rinit(:,:,:,neuronTyping.cellID(:,1,4) == subjecti);
kappas = repmat(permute(rrModel.kappas,[1,5,2,3,4]),[1,size(Y,1),1,1,1]);
input = repmat(permute(rrModel.theoreticalInput,[1,5,2,3,4]),[1,size(Y,1),1,1,1]);
inputsTheory = rrModel.theoreticalInput;
t = rrModel.modelFEF.t;
data_t = neuronTyping.initCoh.neuron_t;

%% Set the vectors
m_r = zeros(size(Y,1),length(rrModel.modelFEF.dimNames));
for di = 1:length(rrModel.modelFEF.dimNames)
    disp(size(m_r))
    switch rrModel.modelFEF.dimNames{di}
        case 'Speed'
            m_r(idx==speedCluster,di) = 1;
        case 'Gain'
            for clusti = 1:neuronTyping.NumClusters
                m_r(idx==clusti,di) = neuronTyping.gainRegression(subjecti).B(clusti);
            end
        otherwise
            error(['Reduced rank model dimension ' rrModel.modelFEF.dimNames{di} ' not recognized!'])
    end
end

I_l = MTtoFEF.vOpt';
for inputi = 1:size(I_l,2)
    I_l_ortho(:,inputi) = I_l(:,inputi) - m_r*inv(m_r'*m_r)*m_r'*I_l(:,inputi);
end

m_r = m_r./vecnorm(m_r);
I_l_ortho = I_l_ortho./vecnorm(I_l_ortho);


%% Estimate R from first principles
Rhat = permute(sum(kappas.*repmat(m_r',[1,1,size(kappas,3),size(kappas,4),size(kappas,5)]),1) + ...
    sum(input(1:2,:,:,:,:).*repmat(I_l_ortho',[1,1,size(input,3),size(input,4),size(input,5)]),1),[3,4,5,2,1]);
Rhat_c = Rhat./max(abs(Rhat),[],[1,2,3]);
R_c = R./max(abs(R),[],[1,2,3]);
R_c = R_c-mean(R_c(data_t<=0,:,:,:),[1,2,3]);
% R_c = (R-mean(R(data_t<=0,:,:,:),'all'))./max(abs(R-mean(R(data_t<=0,:,:,:),'all')),[],[1,2,3]);

%% Estimate R from regression against inputs/dynamics
I_ortho = [I_l zeros(size(I_l_ortho,1),1)];
% angles = zeros(size(inputsTheory,1),length(rrModel.modelFEF.dimNames));
    X = nan(size(inputsTheory,2)*size(inputsTheory,3)*size(inputsTheory,4),size(inputsTheory,1)+size(kappas,1)+1);
    Xin = nan(size(inputsTheory,2)*size(inputsTheory,3)*size(inputsTheory,4),1+size(kappas,1)+1);
    Xout = nan(size(Xin,1),size(R,4));
    Xunexplained = nan(size(Xin,1),size(R,4));
    Xexplained = permute(sum(input(1:2,:,:,:,:).*repmat(I_ortho(:,1:size(I_ortho,2)-1)',[1,1,size(input,3),size(input,4),size(input,5)]),1),[3,2,4,5,1]);
    
    ind = 1;
    for si = 1:size(inputsTheory,3)
        for ci = 1:size(inputsTheory,4)
            X(ind:ind+size(inputsTheory,2)-1,1:size(inputsTheory,1)) = permute(inputsTheory(:,:,si,ci),[2,1,3,4]);
            X(ind:ind+size(inputsTheory,2)-1,size(inputsTheory,1)+1:end-1) = permute(kappas(:,1,:,si,ci),[3,1,2,4,5]);
            X(ind:ind+size(inputsTheory,2)-1,end) = 1;
            Xin(ind:ind+size(inputsTheory,2)-1,1) = permute(inputsTheory(3,:,si,ci),[2,1,3,4]);
            Xin(ind:ind+size(inputsTheory,2)-1,2:end-1) = permute(kappas(:,1,:,si,ci),[3,1,2,4,5]);
            Xin(ind:ind+size(inputsTheory,2)-1,end) = 1;
            Xout(ind:ind+size(inputsTheory,2)-1,:) = squeeze(R_c(data_t<=t(end),si,ci,:));
            Xunexplained(ind:ind+size(inputsTheory,2)-1,:) = Xout(ind:ind+size(inputsTheory,2)-1,:)-Xexplained(:,:,si,ci);
            ind = ind + size(inputsTheory,2);
        end
    end
    
    % Collect vectors
    B = inv(Xin'*Xin + ridgeLambda*eye(size(Xin'*Xin)))*Xin'*Xunexplained;
    B = [I_ortho(:,1:size(I_ortho,2)-1)'; B];
    B2 = (B'./vecnorm(B'))';
    
    % Measure angles between input vectors and
    angles = acosd(B2(1:size(inputsTheory,1),:)*B2(size(inputsTheory,1)+1:end-1,:)');
    
    % Orthogonalize input space to recurrent space
    I_l = B(1:size(inputsTheory,1),:)';
    m_r_temp = B(size(inputsTheory,1)+1:end-1,:)';
    for inputi = 1:size(I_ortho,2)
        I_ortho(:,inputi) = I_l(:,inputi) - m_r_temp*inv(m_r_temp'*m_r_temp)*m_r_temp'*I_l(:,inputi);
    end
    
%     I_ortho = I_ortho./vecnorm(I_ortho);

while any(abs(angles(:)-90) > tol)    
    X = nan(size(inputsTheory,2)*size(inputsTheory,3)*size(inputsTheory,4),size(inputsTheory,1)+size(kappas,1)+1);
    Xin = nan(size(inputsTheory,2)*size(inputsTheory,3)*size(inputsTheory,4),size(kappas,1)+1);
    Xout = nan(size(Xin,1),size(R,4));
    Xunexplained = nan(size(Xin,1),size(R,4));
    Xexplained = permute(sum(input(1:3,:,:,:,:).*repmat(I_ortho',[1,1,size(input,3),size(input,4),size(input,5)]),1),[3,2,4,5,1]);
    
    ind = 1;
    for si = 1:size(inputsTheory,3)
        for ci = 1:size(inputsTheory,4)
            X(ind:ind+size(inputsTheory,2)-1,1:size(inputsTheory,1)) = permute(inputsTheory(:,:,si,ci),[2,1,3,4]);
            X(ind:ind+size(inputsTheory,2)-1,size(inputsTheory,1)+1:end-1) = permute(kappas(:,1,:,si,ci),[3,1,2,4,5]);
            X(ind:ind+size(inputsTheory,2)-1,end) = 1;
            Xin(ind:ind+size(inputsTheory,2)-1,1:end-1) = permute(kappas(:,1,:,si,ci),[3,1,2,4,5]);
            Xin(ind:ind+size(inputsTheory,2)-1,end) = 1;
            Xout(ind:ind+size(inputsTheory,2)-1,:) = squeeze(R_c(data_t<=t(end),si,ci,:));
            Xunexplained(ind:ind+size(inputsTheory,2)-1,:) = Xout(ind:ind+size(inputsTheory,2)-1,:)-Xexplained(:,:,si,ci);
            ind = ind + size(inputsTheory,2);
        end
    end
    
    % Collect vectors
    B = inv(Xin'*Xin + ridgeLambda*eye(size(Xin'*Xin)))*Xin'*Xunexplained;
    B = [I_ortho'; B];
    B2 = (B'./vecnorm(B'))';
    
    % Measure angles between input vectors and
    angles = acosd(B2(1:size(inputsTheory,1),:)*B2(size(inputsTheory,1)+1:end-1,:)');
    
    % Orthogonalize input space to recurrent space
    I_l = B(1:size(inputsTheory,1),:)';
    m_r_temp = B(size(inputsTheory,1)+1:end-1,:)';
    for inputi = 1:size(I_ortho,2)
        I_ortho(:,inputi) = I_l(:,inputi) - m_r_temp*inv(m_r_temp'*m_r_temp)*m_r_temp'*I_l(:,inputi);
    end
    
end


Xhat = X*B;

CC = corrcoef([Xout,Xhat]);
cc = diag(CC(size(Xout,2)+1:end,1:size(Xout,2)));

for clusti = 1:neuronTyping.NumClusters
    mean_cc_by_cluster(clusti) = mean(cc(idx==clusti));
    ste_cc_by_cluster(clusti) = std(cc(idx==clusti))/sqrt(sum(idx==clusti));
end

%% Set covariance matrix
C = cov(B');
Sigma = nan(size(inputsTheory,1)+size(kappas,1)+size(kappas,1));

% Diagonal elements
for i = 1:size(inputsTheory,1)
    Sigma(i,i) = sqrt(rrModel.modelFEF.sigmas(size(kappas,1)+i));     % Input
end
for i = size(inputsTheory,1)+1:size(inputsTheory,1)+size(kappas,1)
    Sigma(i,i) = NaN;       % Selection vectors are unknown
end
ind = 0;
for i = size(inputsTheory,1)+size(kappas,1)+1: size(inputsTheory,1)+size(kappas,1)+size(kappas,1)
    ind = ind+1;
    Sigma(i,i) = sqrt(rrModel.modelFEF.sigmas(ind));
end

% Off diagonal elements
for i = 1:size(inputsTheory,1)
    for j = 2:size(inputsTheory,1)
        if i ~= j
            Sigma(i,j) = C(i,j);        % Input vector covariance estimate
            Sigma(j,i) = C(j,i);
        end
    end
end
indi = size(kappas,1);
for i = 1:size(inputsTheory,1)
    indi = indi+1;
    indj = 0;
    for j = size(inputsTheory,1)+1:size(inputsTheory,1)+size(kappas,1)
        indj = indj+1;
        if i ~= j
            Sigma(i,j) = rrModel.modelFEF.overlaps(indj,indi);     % Input vector and selection vector covariance
            Sigma(j,i) = rrModel.modelFEF.overlaps(indj,indi);     % Input vector and selection vector covariance
        end
    end
end

indi = 0;
for i = size(inputsTheory,1)+1:size(inputsTheory,1)+size(kappas,1)
    indi = indi+1;
    indj = 0;
    for j = size(inputsTheory,1)+size(kappas,1)+1:size(inputsTheory,1)+2*size(kappas,1)
        indj = indj+1;
        if i ~= j
            Sigma(i,j) = rrModel.modelFEF.overlaps(indj,indi);     % output vector and selection vector covariance
            Sigma(j,i) = rrModel.modelFEF.overlaps(indj,indi);     % output vector and selection vector covariance
        end
    end
end

Sigma(size(inputsTheory,1)+2*size(kappas,1)-1,size(inputsTheory,1)+2*size(kappas,1)) = C(size(inputsTheory,1)+size(kappas,1),size(inputsTheory,1)+size(kappas,1)-1);
Sigma(size(inputsTheory,1)+2*size(kappas,1),size(inputsTheory,1)+2*size(kappas,1)-1) = C(size(inputsTheory,1)+size(kappas,1),size(inputsTheory,1)+size(kappas,1)-1);

Sigma(1:size(inputsTheory,1),size(inputsTheory,1)+size(kappas,1)+1:size(inputsTheory,1)+2*size(kappas,1)) = 0;      % Input and output vectors are assumed to be uncorrelated
Sigma(size(inputsTheory,1)+size(kappas,1)+1:size(inputsTheory,1)+2*size(kappas,1),1:size(inputsTheory,1)) = 0;

%% Estimate mean and covariance of each element of the selection vectors
mu = zeros(size(inputsTheory,1)+size(kappas,1)*2,1);
x_indices = false(size(inputsTheory,1)+size(kappas,1)*2,1);
x_indices(1:size(inputsTheory,1)) = true;
x_indices(size(inputsTheory,1)+size(kappas,1)+1:end) = true;

[mu_, Sig_, muUn, SigUnObs, SigObsObs] = conditionalGaussian(mu,Sigma,B(1:end-1,:),'x_indices',x_indices);

%% Functional response topography collapsed onto 1D ring
thetas = atan2(Y(:,2)-mean(Y(:,2)),Y(:,1)-mean(Y(:,1)));
[~,thetaSort] = sort(thetas);

%% Plotting

%% The population vectors represented in the functional topology
figure
subplot(2,2,1)
scatter(Y(m_r(:,1)>0,1),Y(m_r(:,1)>0,2),abs(m_r(m_r(:,1)>0,1))*1000,m_r(m_r(:,1)>0,1),'filled')
axis([-30 30 -30 30])
subplot(2,2,2)
scatter(Y(:,1),Y(:,2),abs(m_r(:,2))*1000,m_r(:,2),'filled')
axis([-30 30 -30 30])
subplot(2,2,3)
scatter(Y(:,1),Y(:,2),abs(I_l_ortho(:,1))*1000,I_l(:,1),'filled')
axis([-30 30 -30 30])
subplot(2,2,4)
scatter(Y(:,1),Y(:,2),abs(I_l_ortho(:,2))*1000,I_l(:,2),'filled')
axis([-30 30 -30 30])

%% The population vectors identified by regression
figure
for regressori = 1:size(B,1)
    subplot(2,size(B,1)+1,regressori)
    scatter(Y(:,1),Y(:,2),abs(B2(regressori,:)')/max(abs(B2),[],'all')*100,B2(regressori,:)','filled')
    caxis([min(B2,[],'all'), max(B2,[],'all')])
    xlabel('tSNE1')
    ylabel('tSNE2')
    axis([-30 30 -30 30])
    
    if regressori < 4
%         A = B2(regressori,:).*input(regressori,:,:,:,:);
        A = B2(regressori,:).*permute(R_c(data_t<=t(end),:,:,:),[5,4,1,2,3]);
    elseif regressori < 6
%         A = B2(regressori,:).*kappas(regressori-3,:,:,:,:);
        A = B2(regressori,:).*permute(R_c(data_t<=t(end),:,:,:),[5,4,1,2,3]);
    else
        A = nan(size(A));
    end
    subplot(2,size(B,1)+1,regressori + size(B,1)+1)
    for si = 1:3
        for ci = 1:3
            plot(t,squeeze(nanmean(A(:,:,:,si,ci),2)),'Color',[speedColors(si,:) cohs(ci)/100],...
                'DisplayName',['speed = ' num2str(speeds(si)) ', coh = ' num2str(cohs(ci))])
            hold on
        end
    end
    xlabel('Time from motion onset (ms)')
    ylabel('Average along this population mode')
end
subplot(2,size(B,1)+1,regressori+1)
scatter(Y(:,1),Y(:,2),abs(cc)*100,cc,'filled')
title('R-values of regression model')
xlabel('tSNE1')
ylabel('tSNE2')
axis([-30 30 -30 30])

subplot(2,size(B,1)+1,regressori+2+size(B,1))
errorbar(1:neuronTyping.NumClusters,mean_cc_by_cluster,ste_cc_by_cluster,'o')
xlabel('Cluster #')
ylabel('Mean correlation coefficient')

%% Predicted firing rates and mean firing rate by cluster
figure
for idxi = 1:neuronTyping.NumClusters
    subplot(3,3,idxi)
    plot(t,squeeze(Rhat_c(:,2,2,idx==idxi)),'Color',[0.6 0.6 0.6])
    hold on
    plot(data_t(data_t<=t(end)),nanmean(R_c(data_t<=t(end),2,2,idx==idxi),4),'k','LineWidth',2)
end

%% Predicted firing rates vs firing rates by cluster
figure
for idxi = 1:neuronTyping.NumClusters
    subplot(3,3,idxi)
    dataTemp = R(data_t<=100,:,:,idx==idxi);
    predictionTemp = Rhat(t<=100,:,:,idx==idxi);
    plot(dataTemp(:),predictionTemp(:),'o')
    hold on
end
    


