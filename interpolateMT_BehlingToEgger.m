function [interpolatedR, Ss, Cs, Ss2, Cs2] = interpolateMT_BehlingToEgger(R,speeds,cohs,speedsFEF,cohsFEF)
%%
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'R')
addRequired(Parser,'speeds')
addRequired(Parser,'cohs')
addRequired(Parser,'speedsFEF')
addRequired(Parser,'cohsFEF')

parse(Parser,R,speeds,cohs,speedsFEF,cohsFEF)

R = Parser.Results.R;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
speedsFEF = Parser.Results.speedsFEF;
cohsFEF = Parser.Results.cohsFEF;

%% Peform interpolation
speedDifferences = any(~ismember(speedsFEF,speeds));
cohDifferences = any(~ismember(cohsFEF,cohs));

[Ss,Cs] = meshgrid(speeds,cohs);
Ss = Ss';
Cs = Cs';
[Ss2,Cs2] = meshgrid(speedsFEF,cohsFEF);
Ss2 = Ss2';
Cs2 = Cs2';

if speedDifferences && cohDifferences
    for neuroni = 1:size(R,4)
        for ti = 1:size(R,1)
            temp = squeeze(R(ti,:,:,neuroni));
            if any(sum(~isnan(temp'),1)>2) && any(sum(~isnan(temp'),2)>2) % Check if there at least two data points on each dimension for interpolation
                % Strip any row or column that doesnt' have at least 2 data points for interpolation
                interpolatedR(ti,:,:,neuroni) = interp2(Ss(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
                    Cs(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',...
                    temp(sum(isnan(temp),2)<3,sum(isnan(temp),1)<3)',Ss2',Cs2','spline')';
            else
                interpolatedR(ti,:,:,neuroni) = nan(size(Ss2));
            end
        end
    end

elseif speedDifferences
    cohInds = ismember(cohsFEF,cohs);
    R = R(:,:,cohInds,:);
    for neuroni = 1:size(R,4)
        for ti = 1:size(R,1)
            for ci = 1:length(cohs)
                interpolatedR(ti,:,ci,neuroni) = interp1(speeds,squeeze(R(ti,:,ci,neuroni)),speedsFEF);
            end
        end
    end
    
elseif cohDifferences
    speedInds = ismember(speedsFEF,speeds);
    R = R(:,speedInds,:,:);
    for neuroni = 1:size(R,4)
        for ti = 1:size(R,1)
            for si = 1:length(speeds)
                interpolatedR(ti,si,:,neuroni) = interp1(cohs,squeeze(R(ti,si,:,neuroni)),cohsFEF);
            end
        end
    end
    
else
    speedInds = ismember(speedsFEF,speeds);
    cohInds = ismember(cohsFEF,cohs);
    warning('No interpolation required, FEF speeds and cohs match MT speeds and cohs.')
    interpolatedR = R(:,speedInds,cohInds,:);
    
end