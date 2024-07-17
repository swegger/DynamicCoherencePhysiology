function [MT, spref, cohsMT, MTnull] = simulateMTdata(mt,t,modelN,varargin)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'mt')
addRequired(Parser,'t')
addRequired(Parser,'modelN')
addParameter(Parser,'speedPrefs',logspace(0,log10(128),1000))
addParameter(Parser,'speedsMT',[2,4,8,16,32])
addParameter(Parser,'cohsMT',[100,30,10])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'normalize',false)
addParameter(Parser,'gaussianApprox',false)
addParameter(Parser,'removeBaseline',false)

parse(Parser,mt,t,modelN,varargin{:})

mt = Parser.Results.mt;
t = Parser.Results.t;
modelN = Parser.Results.modelN;
speedPrefs = Parser.Results.speedPrefs;
speedsMT = Parser.Results.speedsMT;
cohsMT = Parser.Results.cohsMT;
opponentMT= Parser.Results.opponentMT;
normalize = Parser.Results.normalize;
gaussianApprox = Parser.Results.gaussianApprox;
removeBaseline = Parser.Results.removeBaseline;

%% Parameters for each measured neuron
foundation = horzcat(mt.foundation{:});
foundationNull = horzcat(mt.foundationNull{:});
for ni = 1:size(mt.modulationSP,1)
    for ci = 1:length(cohsMT)
        modulationH(:,ni,ci) = mt.modulationH{ni}{ci};
        modulationSP(:,ni,ci) = mt.modulationSP{ni}{ci};
        modulationHnull(:,ni,ci) = mt.modulationHnull{ni}{ci};
        modulationSPnull(:,ni,ci) = mt.modulationSPnull{ni}{ci};
    end
end

if removeBaseline
    foundation(6,:) = 0;
    foundationNull(6,:) = 0;
    modulationH(6,:,:) = 0;
    modulationHnull(6,:,:) = 0;
    modulationSP(6,:,:) = 0;
    modulationSPnull(6,:,:) = 0;
end

if gaussianApprox
   mu.foundation = mean(foundation,2);
   Sig.foundation = cov(foundation');
   mu.foundationNull = mean(foundationNull,2);
   Sig.foundationNull = cov(foundationNull');
   
   for ci = 1:length(cohsMT)
       mu.modulationH(:,ci) = mean(modulationH(:,:,ci),2);
       Sig.modulationH(:,:,ci) = cov(modulationH(:,:,ci)');
       
       mu.modulationSP(:,ci) = mean(modulationSP(:,:,ci),2);
       Sig.modulationSP(:,:,ci) = cov(modulationSP(:,:,ci)');       
       
       mu.modulationHnull(:,ci) = mean(modulationHnull(:,:,ci),2);
       Sig.modulationHnull(:,:,ci) = cov(modulationHnull(:,:,ci)');
       
       mu.modulationSPnull(:,ci) = mean(modulationSPnull(:,:,ci),2);
       Sig.modulationSPnull(:,:,ci) = cov(modulationSPnull(:,:,ci)');
   end
end

%% Sampling
vp = randsample(size(mt.modulationSP,1),modelN,true);
spref = randsample(speedPrefs,modelN,true);

%% Simulate each model neuron
MT = nan(length(t),length(speedsMT),length(cohsMT),modelN);
for ni = 1:modelN
    disp(['File ' num2str(ni) ' of ' num2str(length(mt))])
    
    MTtemp = nan(length(t),length(speedsMT),length(cohsMT));
    
    
    for si = 1:length(speedsMT)
        d = abs(log2(spref(ni)) - log2(speedsMT(si)));
        for ci = 1:length(cohsMT)
%             if opponentMT
                if gaussianApprox
                    foundationTemp = mvnrnd(mu.foundation,Sig.foundation);
                    modulationHtemp = mvnrnd(mu.modulationH(:,ci),Sig.modulationH(:,:,ci));
                    modulationSPtemp = mvnrnd(mu.modulationSP(:,ci),Sig.modulationSP(:,:,ci));
                    
                    foundationNulltemp = mvnrnd(mu.foundationNull,Sig.foundationNull);
                    modulationHnullTemp = mvnrnd(mu.modulationHnull(:,ci),Sig.modulationHnull(:,:,ci));
                    modulationSPnullTemp = mvnrnd(mu.modulationSPnull(:,ci),Sig.modulationSPnull(:,:,ci));
                    
                    MTtempPref(:,si,ci) = modelneuron(t,foundationTemp .* ...
                        modulationHtemp + ...
                        d.*modulationSPtemp);
                    MTtempNull(:,si,ci) = modelneuron(t,foundationNulltemp .* ...
                        modulationHnullTemp + ...
                        d.*modulationSPnullTemp);
                    
                else
                    MTtempPref(:,si,ci) = modelneuron(t,foundation(:,vp(ni)) .* ...
                        modulationH(:,vp(ni),ci) + ...
                        d.*modulationSP(:,vp(ni),ci));
                    MTtempNull(:,si,ci) = modelneuron(t,foundationNull(:,vp(ni)) .* ...
                        modulationHnull(:,vp(ni),ci) + ...
                        d.*modulationSPnull(:,vp(ni),ci));
                end
%             else
%                 if gaussianApprox
%                     foundationTemp = mvnrnd(mu.foundation,Sig.foundation);
%                     modulationHtemp = mvnrnd(mu.modulationH(:,ci),Sig.modulationH(:,:,ci));
%                     modulationSPtemp = mvnrnd(mu.modulationSP(:,ci),Sig.modulationSP(:,:,ci));   
%                     MTtemp(:,si,ci) = modelneuron(t,foundationTemp .* ...
%                         modulationHtemp + ...
%                         d.*modulationSPtemp);      
%                     
%                 else
%                     MTtemp(:,si,ci) = modelneuron(t,foundation(:,vp(ni)) .* ...
%                         modulationH(:,vp(ni),ci) + ...
%                         d.*modulationSP(:,vp(ni),ci));
%                 end
%             end
        end
    end
    if normalize && opponentMT
        MTtemp = (MTtempPref-MTtempNull)/max(MTtempPref,[],[1,2,3]);
        MTtempNull = (MTtempNull-MTtempPref)/max(MTtempPref,[],[1,2,3]);
    elseif normalize
        MTtempNull = MTtempNull/max(MTtemp,[],[1,2,3]);
        MTtemp = MTtemp/max(MTtemp,[],[1,2,3]);
    elseif opponentMT
        MTtemp = MTtempPref-MTtempNull;
        MTtempNull = MTtempNull-MTtempPref;
    else
        MTtemp = MTtempPref;
    end
    
    MT(:,:,:,ni) = MTtemp;
    MTnull(:,:,:,ni) = MTtempNull;
end

%% Functions

%% modelneuron
function out = modelneuron(x,p)
    %%
    out = gauss4(x,p(1),90,40) + gauss4(x,p(2),170,60) + gauss4(x,p(3),300,75) + gauss4(x,p(4),500,150) + gauss4(x,p(5),900,300) + p(6);
    
%% gauss4
function out = gauss4(x,amp,mu,sig)
    %%
    out = amp.*exp( -((x-mu)./sig).^2 );