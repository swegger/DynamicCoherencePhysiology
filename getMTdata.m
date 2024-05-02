function [MT, spref, swidth, dirMT] = getMTdata(mt,varargin)
%%
%
%
%
%%

%% Defaults
speedPrefOpts_default.tWin = [40,120];
speedPrefOpts_default.P0 = [16,1];
speedPrefOpts_default.ub = [254,128];
speedPrefOpts_default.lb = [0, 0];
speedPrefOpts_default.c = NaN;
speedPrefOpts_default.s = NaN;
speedPrefOpts_default.d = 0;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'mt')
addParameter(Parser,'speedsMT',[2,4,8,16,32])
addParameter(Parser,'cohsMT',[10 30 70 100])
addParameter(Parser,'directionsMT',0)
addParameter(Parser,'opponentMT',false)
addParameter(Parser,'normalize',false)
addParameter(Parser,'sprefFromFit',true)
addParameter(Parser,'speedPrefOpts',speedPrefOpts_default)
addParameter(Parser,'checkMTFit',false)

parse(Parser,mt,varargin{:})

mt = Parser.Results.mt;
speedsMT = Parser.Results.speedsMT;
cohsMT = Parser.Results.cohsMT;
directionsMT = Parser.Results.directionsMT;
opponentMT= Parser.Results.opponentMT;
normalize = Parser.Results.normalize;
sprefFromFit = Parser.Results.sprefFromFit;
speedPrefOpts = Parser.Results.speedPrefOpts;
checkMTFit = Parser.Results.checkMTFit;


%% Get mean of each unit
MT = [];
spref = [];
swidth = [];
dirMT = [];
if opponentMT
    dN = 1;
else
    dN = length(directionsMT);
end
for di = 1:dN
    for filei = 1:length(mt)
        disp(['File ' num2str(filei) ' of ' num2str(length(mt))])
    
        MTtemp = nan(length(mt{filei}.neuron_t),length(speedsMT),length(cohsMT));
        
        
        for si = 1:length(speedsMT)
            for ci = 1:length(cohsMT)
                if opponentMT
                    [~,condLogical] = trialSort(mt{filei},directionsMT(1),speedsMT(si),NaN,cohsMT(ci));
                    [~,condLogicalNull] = trialSort(mt{filei},directionsMT(2),speedsMT(si),NaN,cohsMT(ci));
                    MTtempPref(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
                    MTtempNull(:,si,ci) = mean(mt{filei}.r(:,condLogicalNull),2);
                else
                    [~,condLogical] = trialSort(mt{filei},directionsMT(di),speedsMT(si),NaN,cohsMT(ci));
                    MTtemp(:,si,ci) = mean(mt{filei}.r(:,condLogical),2);
                end
            end
        end
        if normalize && opponentMT
            MTtemp = (MTtempPref-MTtempNull)/max(MTtempPref,[],[1,2,3]);
        elseif normalize
            MTtemp = MTtemp/max(MTtemp,[],[1,2,3]);
        elseif opponentMT
            MTtemp = MTtempPref-MTtempNull;
        end
        
        if sprefFromFit
            [mu,sig,~,normR,~] = fitSpeedTuning(mt{filei},'P0',speedPrefOpts.P0,...
                'ub',speedPrefOpts.ub,'lb',speedPrefOpts.lb,...
                'c',speedPrefOpts.c,'s',speedPrefOpts.s,'d',speedPrefOpts.d,...
                'tWin',speedPrefOpts.tWin);
            spref = [spref mu];
            swidth = [swidth sig];
            
            if checkMTFit
                s = linspace(min(speedsMT)-0.1*min(speedsMT),max(speedsMT)*1.1,20);
                h = figure;
                semilogx(mt{filei}.speeds,normR,'o')
                hold on
                semilogx(s,speedTuning(mt{filei},s,mu,sig))
                ax = axis;
                text(ax(1),0.95*ax(4),['\mu = ' num2str(mu) ', \sig = ' num2str(sig)])
                input('Press enter to continue ')
                close(h);
            end
        else
            spref = [spref speedsMT(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100)) == ...
                max(nansum(MTtemp(mt{filei}.neuron_t >= 50 & mt{filei}.neuron_t <= 150,:,cohs == 100))))];
            swidth = [swidth, nan(size(MTtemp,2))];
        end
        dirMT = [dirMT, directionsMT(di)];
        MT = cat(4,MT,MTtemp);
    end
end