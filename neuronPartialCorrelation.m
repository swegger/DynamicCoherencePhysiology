function [RHO, PVAL, Z] = neuronPartialCorrelation(dcp,varargin)
%% neuronPartialCorrelation
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dcp')
addParameter(Parser,'sourceDirectory','/mnt/Lisberger/Experiments/DynamicCoherencePhysiology/data/Aristotle')
addParameter(Parser,'regressionModels',{'speed','coh','gain'})
addParameter(Parser,'speeds',[5; 10; 20])
addParameter(Parser,'cohs',[20; 60; 100])
addParameter(Parser,'gains',[])
addParameter(Parser,'directions',[0 180])
addParameter(Parser,'tWin',[750 750])
addParameter(Parser,'tGainMeasurement',750)
addParameter(Parser,'gainFitInd',1)
addParameter(Parser,'saveFigures',false)
addParameter(Parser,'saveResults',false)

parse(Parser,dcp,varargin{:})

dcp = Parser.Results.dcp;
sourceDirectory = Parser.Results.sourceDirectory;
regressionModels = Parser.Results.regressionModels;
speeds = Parser.Results.speeds;
cohs = Parser.Results.cohs;
gains = Parser.Results.gains;
directions = Parser.Results.directions;
tWin = Parser.Results.tWin;
tGainMeasurement = Parser.Results.tGainMeasurement;
gainFitInd = Parser.Results.gainFitInd;
saveFigures = Parser.Results.saveFigures;
saveResults = Parser.Results.saveResults;

%% Preliminary
[Cohs, Spds] = meshgrid(cohs,speeds);

% Colors
speedColors = projectColorMaps_coh('speeds','sampleDepth',length(speeds),'sampleN',length(speeds));
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

%% Remove data that doesn't pass cutoff
Rinit = Rinit(:,:,:,passCutoff);
cellID = cellID(passCutoff,:,:);

%% Fit linear models to each neuron
for uniti = 1:size(Rinit,4)
    for modeli = 1:length(regressionModels)
        
        switch regressionModels{modeli}
            
            case 'speed'
                % Get firing rate data
                r = Rinit(taccept,:,:,uniti);
                r = reshape(permute(r,[2,3,1]),[numel(Spds)*sum(taccept),1]);
                mfit(modeli,uniti).R = r;
                
                % Build regressor matrix
                mfit(modeli,uniti).X = repmat([Spds(:),ones(numel(Spds),1)],[sum(taccept),1]);
                df(modeli) = numel(Spds)*sum(taccept)-1-1;
                
            case 'coh'
                % Get firing rate data
                r = Rinit(taccept,:,:,uniti);
                r = reshape(permute(r,[2,3,1]),[numel(Cohs)*sum(taccept),1]);
                mfit(modeli,uniti).R = r;
                
                % Build regressor matrix
                mfit(modeli,uniti).X = repmat([Cohs(:),ones(numel(Cohs),1)],[sum(taccept),1]);
                df(modeli) = numel(Cohs)*sum(taccept)-1-1;
                
            case 'gain'
                if isempty(gains)
                    error('No gain data provided!')
                else
                    gtemp = gains(:,:,gainFitInd);
                    
                    % Get firing rate data
                    r = Rinit(taccept,:,:,uniti);
                    r = reshape(permute(r,[2,3,1]),[numel(gtemp)*sum(taccept),1]);
                    mfit(modeli,uniti).R = r;
                    
                    % Build regressor matrix
                    mfit(modeli,uniti).X = repmat([gtemp(:),ones(numel(gtemp),1)],[sum(taccept),1]);
                    df(modeli) = numel(Spds)*sum(taccept)-1-1;
                end
                
            case 'speed&coh'
                % Get firing rate data
                r = Rinit(taccept,:,:,uniti);
                r = reshape(permute(r,[2,3,1]),[numel(Spds)*sum(taccept),1]);
                mfit(modeli,uniti).R = r;
                
                % Build regressor matrix
                mfit(modeli,uniti).X = repmat([Spds(:),Cohs(:),ones(numel(Spds),1)],[sum(taccept),1]);
                df(modeli) = numel(Spds)*sum(taccept)-2-1;
                
            case 'speed&gain'
                if isempty(gains)
                    error('No gain data provided!')
                else
                    % Get firing rate data
                    r = Rinit(taccept,:,:,uniti);
                    r = reshape(permute(r,[2,3,1]),[numel(Spds)*sum(taccept),1]);
                    mfit(modeli,uniti).R = r;
                    
                    % Build regressor matrix
                    mfit(modeli,uniti).X = repmat([Spds(:),gains(:),ones(numel(Spds),1)],[sum(taccept),1]);
                    df(modeli) = numel(Spds)*sum(taccept)-2-1;
                end
                
            case 'coh&gain'
                
                if isempty(gains)
                    error('No gain data provided!')
                else
                    % Get firing rate data
                    r = Rinit(taccept,:,:,uniti);
                    r = reshape(permute(r,[2,3,1]),[numel(Spds)*sum(taccept),1]);
                    mfit(modeli,uniti).R = r;
                    
                    % Build regressor matrix
                    mfit(modeli,uniti).X = repmat([Cohs(:),gains(:),ones(numel(Cohs),1)],[sum(taccept),1]);
                    df(modeli) = numel(Cohs)*sum(taccept)-2-1;
                end
                
            case 'speed&coh&gain'
                
                if isempty(gains)
                    error('No gain data provided!')
                else
                    % Get firing rate data
                    r = Rinit(taccept,:,:,uniti);
                    r = reshape(permute(r,[2,3,1]),[numel(Spds)*sum(taccept),1]);
                    mfit(modeli,uniti).R = r;
                    
                    % Build regressor matrix
                    mfit(modeli,uniti).X = repmat([Spds(:),Cohs(:),gains(:),ones(numel(Spds),1)],[sum(taccept),1]);
                    df(modeli) = numel(Spds)*sum(taccept)-3-1;
                end
                
            otherwise
                error('Regression model not recognized!')
        end
        
        % Fit model
        [mfit(modeli,uniti).B, mfit(modeli,uniti).BINT, mfit(modeli,uniti).RVAL, mfit(modeli,uniti).RINT, mfit(modeli,uniti).STATS] = ...
            regress(mfit(modeli,uniti).R,mfit(modeli,uniti).X);
        
        % Predict firing rates
        mresult(modeli,uniti).R(:,1) = mfit(modeli,uniti).B'*mfit(modeli,uniti).X';
        
    end
    
    % Perform partial correlation analysis of predicted rates and actual rates
    Rresults = [mfit(1,uniti).R horzcat(mresult(:,uniti).R)];
    [rho, pval] = partialcorr(Rresults);
    
    RHO(:,uniti) = rho(2:end,1);
    PVAL(:,uniti) = pval(2:end,1);
    
    % Perform Fisher's r-to-Z transformation
    for compi = 1:length(regressionModels)
        Z(compi,uniti) = 0.5*log( (1+RHO(compi,uniti))/(1-RHO(compi,uniti)) ) / sqrt(1/df(compi));
    end

end

%% Find significantly different cells
for modeli = 1:length(regressionModels)
    for modelj = 1:length(regressionModels)
        mFavored{modeli,modelj} = find(Z(modeli,:)-Z(modelj,:) > 1.28 & Z(modeli,:) > 1.28);
    end
end

%% Plot results
model1_vs_model2 = figure;
plot(Z(1,:),Z(2,:),'o')
hold on
axis equal
uh = plotUnity;
axis square
ax = axis;
lineProps.LineStyle = '-';
lineProps.Color = 'k';
plotVertical(1.28,'MinMax',[ax(3) 0],'lineProperties',lineProps);
plotHorizontal(1.28,'MinMax',[ax(1) 0],'lineProperties',lineProps);
plot(linspace(0,ax(2),10),linspace(0,ax(2),10)+1.28,'k-')
plot(linspace(1.28,ax(2),10),linspace(1.28,ax(2),10)-1.28,'k-')
axis(ax)
uh.Visible = 'off';
xlabel(['Z-transformed component correlation for ' regressionModels{1}])
ylabel(['Z-transformed component correlation for ' regressionModels{2}])

%% Plot model 2 favored over model 1
for uniti = 1:length(mFavored{2,1})
    gainFavored(uniti) = figure('Name',['Unit ' num2str(cellID(mFavored{2,1}(uniti),1,1)) ' ' num2str(cellID(mFavored{2,1}(uniti),1,2))],'Position', [445 100+uniti*2 1784 420]);
    for si = 1:length(speeds)
        axh(uniti,si) = subplot(1,length(speeds),si);
        for ci = 1:length(cohs)
            plot(initCoh.neuron_t,1000*Rinit(:,si,ci,mFavored{2,1}(uniti)),'Color',(100-cohs(ci))/100*ones(1,3))
            hold on
        end
        axis tight
        plotVertical(tWin);
        ylabel('Spikes/s')
        xlabel('Time from motion onset (ms)')
        title(['Speed = ' num2str(speeds(si))])
    end
    linkaxes(axh(uniti,:))
end

%% Plot average of neurons that favor model 2 over model 1
gradAveModel1_vs_Model2 = figure('Name','Grand average of model 2 favored cells','Position', [445 100+uniti*2 1784 420]);
grandAve = nanmean(1000*Rinit,4);
grandAveFavored = nanmean(1000*Rinit(:,:,:,mFavored{2,1}),4);
for si = 1:length(speeds)
    axh2(si) = subplot(1,length(speeds),si);
    for ci = 1:length(cohs)
        plot(initCoh.neuron_t,grandAveFavored(:,si,ci),'Color',(100-cohs(ci))/100*ones(1,3))
        hold on
    end
    axis tight
    ylabel('Spikes/s')
    xlabel('Time from motion onset (ms)')
    title(['Speed = ' num2str(speeds(si))])
end
linkaxes(axh2)
for si = 1:length(speeds)
    subplot(axh2(si))
    plotVertical(tGainMeasurement);
end

grandAveGainRepresenation = figure('Name',...
    ['Representation of gain by average of' regressionModels{2} ' preferring neurons'],...
    'Position',[447 614 1598 420]);
for gi = 1:length(tGainMeasurement)
    subplot(1,length(tGainMeasurement)+1,gi)
    for si = 1:length(speeds)
        plot(gains(si,:,gi),squeeze(grandAveFavored(initCoh.neuron_t == tGainMeasurement(gi),si,:)),...
            'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si)) ', gain favoring'])
        hold on
        plot(gains(si,:,gi),squeeze(grandAve(initCoh.neuron_t == tGainMeasurement(gi),si,:)),...
            'o-','Color',speedColors(si,:),...
            'DisplayName',['Speed = ' num2str(speeds(si)) ', all neurons'])
    end
    xlabel('Behavioral gain (unitless)')
    ylabel('Mean spikes/s')
    title(['For gain at t = ' num2str(tGainMeasurement(gi))])
end

subplot(1,length(tGainMeasurement)+1,length(tGainMeasurement)+1)
for si = 1:length(speeds)
    plot(gains(si,:,gainFitInd),squeeze(nanmean(grandAveFavored(taccept,si,:),1)),...
        'o-','Color',speedColors(si,:),'MarkerFaceColor',speedColors(si,:),...
        'DisplayName',['Speed = ' num2str(speeds(si)) ', gain favoring'])
    hold on
    plot(gains(si,:,gainFitInd),squeeze(nanmean(grandAve(taccept,si,:),1)),...
        'o-','Color',speedColors(si,:),...
        'DisplayName',['Speed = ' num2str(speeds(si)) ', all neurons'])
end
xlabel('Behavioral gain (unitless)')
ylabel('Mean spikes/s')
title(['For gain at t = ' num2str(tWin)])

%% Save figures
if saveFigures
    
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/mat/' dcp{1}.sname ...
        '/partialCorrelation/' datestr(now,'yyyymmdd')];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    
    savefig(model1_vs_model2,[saveLocation '/zscores_' regressionModels{1} '_vs_' regressionModels{2} '.fig'])
    savefig(gradAveModel1_vs_Model2 ,[saveLocation '/averagePSTH_grainPreferring.fig'])
    savefig(grandAveGainRepresenation ,[saveLocation '/gainRepresenationByGainPreferring.fig'])
    
    for uniti = 1:length(mFavored{2,1})
        savefig(gainFavored(uniti),[saveLocation '/' regressionModels{2} 'PreferringUnit_' num2str(cellID(mFavored{2,1}(uniti),1,1)) '_' num2str(cellID(mFavored{2,1}(uniti),1,2)) '.fig'])
    end
end

%% Saving
if saveResults
    saveLocation = ['/mnt/Lisberger/Manuscripts/FEFphysiology/Subprojects/FEFinitiation/' dcp{1}.sname ...
        '/partialCorrelationResults'];
    if ~exist(saveLocation,'dir')
        mkdir(saveLocation)
    end
    save([saveLocation '/partialCorrelationResults' datestr(now,'yyyymmdd')],'-v7.3')
end