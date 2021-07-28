function plotInitSeqComp(initCoh,dynamicCoh,unit,windows,sensitivity,m,ste,varargin)
%% plotInitSeqComp
%
%
%%

%% Defaults
saveOpts_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'initCoh')
addRequired(Parser,'dynamicCoh')
addRequired(Parser,'unit')
addRequired(Parser,'windows')
addRequired(Parser,'sensitivity')
addRequired(Parser,'m')
addRequired(Parser,'ste')
addParameter(Parser,'t',-100:1600)
addParameter(Parser,'saveOpts',saveOpts_default)
addParameter(Parser,'closeFig',false)

parse(Parser,initCoh,dynamicCoh,unit,windows,sensitivity,m,ste,varargin{:})

initCoh = Parser.Results.initCoh;
dynamicCoh = Parser.Results.dynamicCoh;
unit = Parser.Results.unit;
windows = Parser.Results.windows;
sensitivity = Parser.Results.sensitivity;
m = Parser.Results.m;
ste = Parser.Results.ste;
t = Parser.Results.t;
saveOpts = Parser.Results.saveOpts;
closeFig = Parser.Results.closeFig;

%% Set variables
unitInd = find(initCoh.unitIndex == unit);
winN = size(windows,1);

%% Plot 
fh = figure('Name',['Unit ' initCoh.datapath(end-8:end) ' # ' num2str(unit)],...
    'Position', [294 121 1900 1174]);
ah(1) = subplot(2,2*winN,1:winN);
plot(t,squeeze(initCoh.R(:,:,:,unitInd))*1000)
axis tight
hold on
for wini = 1:winN
    errorbar(mean(windows(wini,:))*ones(winN,1),m.init(:,wini)*1000/100,ste.init(:,wini)*1000/100,'ko')
end
xlabel('Time from motion onset')
ylabel('sp/s')
ax(1,:) = axis;

cohs = unique(initCoh.coh);
cs = linspace(min(cohs)-10,max(cohs)+10,100);
for wini = 1:winN
    ah2(wini) = subplot(2,2*winN,2*winN+wini);
    errorbar(cohs,m.init(:,wini)*1000/100,ste.init(:,wini)*1000/100,'ko')
    hold on
    plot(cs,1000/100*exp(sensitivity.init.theta(1,wini)*cs + sensitivity.init.theta(2,wini)),'k-')
    xlabel('Coherence')
    ylabel('sp/s')
    title([num2str(windows(wini,1)) ' to ' num2str(windows(wini,2)) ' ms'])
    ax2(wini,:) = axis;
end
legend({'Data','LNP'})

ah(2) = subplot(2,2*winN,winN+1:winN*2);
if unitInd > size(dynamicCoh.R,3)
    warning(['dynamicCoh.R for unit ' num2str(dynamicCoh.unitIndex(unitInd)) ' does not exist...'])
else
    plot(t,squeeze(dynamicCoh.R(:,:,unitInd))*1000)
end
axis tight
hold on
for wini = 1:winN
    errorbar(mean(windows(wini,:))*ones(5,1),m.dyn(:,wini)*1000/100,ste.dyn(:,wini)*1000/100,'ko')
end
xlabel('Time from motion onset')
ylabel('sp/s')
ax(2,:) = axis;
axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
for wini = 1:winN
    plotVertical(windows(wini,:));
end
legend({'60->100->20','60->20->100','60->100','60->20','60'})
% mymakeaxis(gca);
subplot(ah(1))
axis([min(ax(:,1)) max(ax(:,2)) min(ax(:,3)) max(ax(:,4))])
for wini = 1:winN
    plotVertical(windows(wini,:));
end
legend({'20%','60%','100%'})
% mymakeaxis(gca);

seqs = unique(dynamicCoh.sequences);
for seqi = 1:length(seqs)
    for wini = 1:winN
        cohmap(seqi,wini) = dynamicCoh.coh(round(mean(windows(wini,:)))+101,seqi);
    end
end
cs = linspace(min(cohs)-10,max(cohs)+10,100);
for wini = 1:winN
    ah2(wini+winN) = subplot(2,2*winN,3*winN+wini);
    errorbar(cohmap(:,wini),m.dyn(:,wini)*1000/100,ste.dyn(:,wini)*1000/100,'ko')
    hold on
    if wini > 1
        plot(cs,1000/100*exp(sensitivity.dyn.theta(1,wini)*cs + sensitivity.dyn.theta(2,wini)),'k-')
    end
    xlabel('Coherence')
    ylabel('sp/s')
    title([num2str(windows(wini,1)) ' to ' num2str(windows(wini,2)) ' ms'])
    ax2(wini+winN,:) = axis;
%     mymakeaxis(gca);
end

for ai = 1:6
    subplot(ah2(ai))
    axis([min(ax2(:,1)) max(ax2(:,2)) min(ax2(:,3)) max(ax2(:,4))])
%     mymakeaxis(gca);
end
legend({'Data','LNP'})

%% Save
if saveOpts.On
    savefig(fh,saveOpts.filename)
end

%% Close
if closeFig
    close(fh)
end