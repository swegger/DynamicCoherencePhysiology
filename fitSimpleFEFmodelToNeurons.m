function modelFEF = fitSimpleFEFmodelToNeurons(dcp,varargin)
%% neuronPartialCorrelation
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'objectFile','allMT_20230419.mat')
addParameter(Parser,'speeds',[2,4,8,16,32])
addParameter(Parser,'cohs',[10 30 70 100])
addParameter(Parser,'speedsFEF',[5,10,20])
addParameter(Parser,'cohsFEF',[20, 60, 100])
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'P0',NaN)
addParameter(Parser,'ub',NaN)
addParameter(Parser,'lb',NaN)
addParameter(Parser,'tWin',[0 900])
addParameter(Parser,'tau',20)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
objectFile = Parser.Results.objectFile;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;
directions = Parser.Results.directions;
P0 = Parser.Results.P0;
ub = Parser.Results.ub;
lb = Parser.Results.lb;
tWin = Parser.Results.tWin;
tau = Parser.Results.tau;

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

%% Get MT data and organize as an input

% Load mtObjs
load(objectFile)

% Get mean of each unit
MT = [];
spref = [];
for filei = 1:length(mt)
    disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
    
    MTtemp = nan(length(mt{filei}.neuron_t),length(speeds),length(cohs));
    
    
    for si = 1:length(speeds)
        for ci = 1:length(cohs)
            [~,condLogical] = trialSort(mt{filei},0,speeds(si),NaN,cohs(ci));
            MTtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
        end
    end
    %if ~any(isnan(Rtemp(:)))
        spref = [spref speeds(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
            max(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100))))];
        MT = cat(4,MT,MTtemp);
    %end
end

mtNeuron_t = mt{filei}.neuron_t;


% Interpolate between coherence speed values from Behling experiments to the speed and coherence values used in Egger experiments
[interpolatedR, ~, ~, ~, ~] = interpolateMT_BehlingToEgger(MT,speeds,cohs,speedsFEF,cohsFEF);
% [Ss,Cs] = meshgrid(speeds,cohs);
% Ss = Ss';
% Cs = Cs';
% [Ss2,Cs2] = meshgrid(speedsFEF,cohsFEF);
% Ss2 = Ss2';
% Cs2 = Cs2';
% for neuroni = 1:size(MT,4)
%     for ti = 1:size(MT,1)
%         temp = squeeze(MT(ti,:,:,neuroni));
%         if any(sum(~isnan(temp'),1)>2) && any(sum(~isnan(temp'),2)>2) % Check if there at least two data points on each dimension for interpolation
%             % Strip any row or column that doesnt' have at least 2 data points for interpolation
%             interpolatedR(ti,:,:,neuroni) = interp2(Ss(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
%                 Cs(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
%                 temp(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',Ss2',Cs2','spline')';
%         else
%              interpolatedR(ti,:,:,neuroni) = nan(size(Ss2));
%         end
%     end
% end

inputs = permute(interpolatedR,[4,1,2,3]);


%% Iterate across neurons and fit data

% Prepare data for fitting
[~,iFEF,iMT] = intersect(fef_t,mtNeuron_t);
RtoFit = Rinit(iFEF,:,:,:);
inputs = inputs(:,iMT,:,:);

if any(isnan(P0))
    P0 = [0.5,10,randn(1,size(inputs,1))];
end
if any(isnan(lb))
    lb = [0, 0, -Inf*ones(1,size(inputs,1))];
end
if any(isnan(ub))
    ub = [1, 100, Inf*ones(1,size(inputs,1))];
end
for neuroni = 1:size(RtoFit,4)
    Rtemp = RtoFit(:,:,:,neuroni);
    R0 = nanmean(Rtemp(1,:,:),[2,3]);
    minimizer = @(P)minimizant(P,R0,tau,inputs,Rtemp);
    [P, fval, exitflag, output, lambda, grad, hessian] = fmincon(minimizer,P0,[],[],[],[],lb,ub);
    
    modelFEF.R0 = R0;
    modelFEF.leakRate(neuroni) = P(1);
    modelFEF.baseLine(neuroni) = P(2);
    modelFEF.W(:,neuroni) = P(3:end);
    modelFEF.fval(:,neuroni) = fval;
end
    
%% Functions

%% Function to fit
function R = SimpleFEFmodel(leakRate,W,baseLine,R0,tau,inputs)
    sz = size(inputs);
    R = nan(sz(2:4));
    R(1,:,:) = repmat(R0,[1,sz(3),sz(4)]);
    for ti = 2:sz(2)
        for si = 1:sz(3)
            for ci = 1:sz(4)
                dR(ti,si,ci) = -leakRate.*(R(ti-1,si,ci)-baseLine) + W*inputs(:,ti,si,ci);
                R(ti,si,ci) = R(ti-1,si,ci) + dR(ti,si,ci)/tau;
            end
        end
    end
    
 %% Objective
function out = minimizant(P,R0,tau,inputs,R)
    leakRate = P(1);
    baseLine = P(2);
    W = P(3:end);
    Rest = SimpleFEFmodel(leakRate,W,baseLine,R0,tau,inputs);
    out = sum((R(:) - Rest(:)).^2);