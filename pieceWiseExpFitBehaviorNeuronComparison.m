params.behavior = nan(length(dynamicCoh),5,4,4);
params.neuron = nan(length(dynamicCoh),5,4,4,10);
response.behavior = nan(length(t),length(dynamicCoh),5,4);
response.neuron = nan(length(t),length(dynamicCoh),5,4,10);
epochN = 5;
seqN = 5;
for sitei = 1:length(dynamicCoh)
    for epochi = 1:epochN
        if epochi < 3
            tempCont = 1;
        else
            tempCont = 6;
        end
        for seqi = 1:seqN
            if ~isempty(dynamicCoh{sitei}.Epoch)
                params.behavior(sitei,seqi,epochi,:) = dynamicCoh{sitei}.Epoch.p(seqi,epochi,tempCont,1,1).eye;
                response.behavior(:,sitei,seqi,epochi) = pieceWiseExpRise(t,dynamicCoh{sitei}.Epoch.p(seqi,epochi,tempCont,1,1).eye);
                for uniti = 1:size(dynamicCoh{sitei}.Epoch.p,5)
                    params.neuron(sitei,seqi,epochi,:,uniti) = dynamicCoh{sitei}.Epoch.p(seqi,epochi,tempCont,1,uniti).neurons;
                    response.neuron(:,sitei,seqi,epochi,uniti) = pieceWiseExpRise(t,dynamicCoh{sitei}.Epoch.p(seqi,epochi,tempCont,1,uniti).neurons);
                    responseReal.neuron(:,sitei,seqi,epochi,uniti) = dynamicCoh{sitei}.Epoch.neurons(:,seqi,epochi,tempCont,1,uniti);
                end
            end
        end
    end
end

%% Plot
unitMap = [];
for sitei = 1:length(dynamicCoh)
    if ~isempty(dynamicCoh{sitei}.Epoch)
        for uniti = 1:size(dynamicCoh{sitei}.Epoch.p,5)
            if ~isnan(response.neuron(1,sitei,1,1,uniti))
                unitMap = [unitMap; sitei uniti];
            end
        end
    end
end
for uniti = 1:size(unitMap,1)
    figure(uniti)
    colors = colormap('lines');
    for epochi = 1:epochN
        subplot(3,epochN,epochi)
        for seqi = 1:seqN
            plot(t,nanmean(response.behavior(:,:,seqi,epochi),2),'-','Color',colors(seqi,:))
            hold on
        end
                
        subplot(3,epochN,epochN+epochi)
        for seqi = 1:seqN
            
            plot(t,1000*nanmean(responseReal.neuron(:,unitMap(uniti,1),seqi,epochi,unitMap(uniti,2)),2),'-','Color',colors(seqi,:))
            hold on
        end
        
        subplot(3,epochN,2*epochN+epochi)
        for seqi = 1:seqN
            
            plot(t,1000*nanmean(response.neuron(:,unitMap(uniti,1),seqi,epochi,unitMap(uniti,2)),2),'-','Color',colors(seqi,:))
            hold on
        end
    end
    title(num2str(unitMap(uniti,:)))
end

%%
figure
goodUnits = [13,2;25,1;26,1;26,1;26,2;29,1;32,2;34,2;36,3];
colors = colormap('lines');
for epochi = 1:epochN
    subplot(2,epochN,epochi)
    for seqi = 1:seqN
        plot(t,nanmean(response.behavior(:,:,seqi,epochi),2),'-','Color',colors(seqi,:))
        hold on
    end
    
    subplot(2,epochN,epochN+epochi)
    for seqi = 1:seqN
        for uniti = 1:size(unitMap,1)
            if any(ismember(unitMap(uniti,:),goodUnits,'rows'))
                AMP = 1000*(params.neuron(unitMap(uniti,1),seqi,epochi,2,unitMap(uniti,2)) - params.neuron(unitMap(uniti,1),seqi,epochi,1,unitMap(uniti,2)));
                plot(t,1000*nanmean(response.neuron(:,unitMap(uniti,1),seqi,epochi,unitMap(uniti,2)),2)/abs(AMP),'-','Color',colors(seqi,:))
            end
            hold on
        end
    end    
end


figure('Name','R_f - R_0')
for epochi = 1:epochN
    subplot(2,epochN,epochi)
    for seqi = 1:seqN
        temp = params.behavior(:,seqi,epochi,2)-params.behavior(:,seqi,epochi,1);
        temp = temp(abs(temp) < 200);
        tempStd = nanstd(temp);
        errorbar(seqi,nanmean(temp),tempStd/sqrt(size(params.behavior,1)),'o','Color',colors(seqi,:))
        hold on
    end
    plotHorizontal(0);
    
    subplot(2,epochN,epochN+epochi)
    for seqi = 1:seqN
        tempparam = [];
        for uniti = 1:size(unitMap,1)
            if any(ismember(unitMap(uniti,:),goodUnits,'rows'))
                tempparam = [tempparam; 1000*(params.neuron(unitMap(uniti,1),seqi,epochi,2,unitMap(uniti,2)) - params.neuron(unitMap(uniti,1),seqi,epochi,1,unitMap(uniti,2)))];
            end
        end
        temp = nanmean(tempparam(abs(tempparam) < 500));
        tempStd = nanstd(tempparam(abs(tempparam) < 500));
        errorbar(seqi,temp,tempStd/sqrt(sum(abs(tempparam) < 500)),'o','Color',colors(seqi,:))
        hold on
%         tempparam2(:,seqi) = tempparam;
    end
    plotHorizontal(0);
%     tempparam2(abs(tempparam2) > 500) = NaN;
%     plot(1:seqN,tempparam2,'.-','Color',colors(seqi,:))
end

figure('Name','\tau')
for epochi = 1:epochN
    subplot(2,epochN,epochi)
    for seqi = 1:seqN
        temp = nanmean(params.behavior(params.behavior(:,seqi,epochi,3) < 500,seqi,epochi,3));
        tempStd = nanstd(params.behavior(params.behavior(:,seqi,epochi,3) < 500,seqi,epochi,3));
        errorbar(seqi,temp,tempStd/sqrt(sum(params.behavior(:,seqi,epochi,3) < 500)),'o','Color',colors(seqi,:))
        hold on
    end
    
    subplot(2,epochN,epochN+epochi)
    for seqi = 1:seqN
        tempparam = [];
        for uniti = 1:size(unitMap,1)
            if any(ismember(unitMap(uniti,:),goodUnits,'rows'))
                tempparam = [tempparam; params.neuron(unitMap(uniti,1),seqi,epochi,3,unitMap(uniti,2))];
            end
        end
        temp = nanmean(tempparam(tempparam < 500));
        tempStd = nanstd(tempparam(tempparam < 500));
        errorbar(seqi,temp,tempStd/sqrt(sum(tempparam < 500)),'o','Color',colors(seqi,:))
        hold on
    end
end

figure('Name','t_0')
for epochi = 1:epochN
    subplot(2,epochN,epochi)
    for seqi = 1:seqN
        temp = nanmean(params.behavior(params.behavior(:,seqi,epochi,4) < 500,seqi,epochi,4));
        tempStd = nanstd(params.behavior(params.behavior(:,seqi,epochi,4) < 500,seqi,epochi,4));
        errorbar(seqi,temp,tempStd/sqrt(sum(params.behavior(:,seqi,epochi,4) < 500)),'o','Color',colors(seqi,:))
        hold on
    end
    
    subplot(2,epochN,epochN+epochi)
    for seqi = 1:seqN
        tempparam = [];
        for uniti = 1:size(unitMap,1)
            if any(ismember(unitMap(uniti,:),goodUnits,'rows'))
                tempparam = [tempparam; params.neuron(unitMap(uniti,1),seqi,epochi,3,unitMap(uniti,2))];
            end
        end
        temp = nanmean(tempparam(tempparam < 500));
        tempStd = nanstd(tempparam(tempparam < 500));
        errorbar(seqi,temp,tempStd/sqrt(sum(tempparam<500)),'o','Color',colors(seqi,:))
        hold on
    end
end

%%
figure
ind = 0;
for parami = 3:4
    ind = ind+1;
    for epochi = 1:epochN
        subplot(2,epochN,(ind-1)*epochN +epochi)
        for seqi = 1:seqN
            for uniti = 1:size(unitMap,1)
                if any(ismember(unitMap(uniti,:),goodUnits,'rows'))
                    plot(params.behavior(unitMap(uniti,1),seqi,epochi,parami),...
                        params.neuron(unitMap(uniti,1),seqi,epochi,parami,unitMap(uniti,2)),...
                        's','Color',colors(seqi,:))
                    hold on
                end
            end
        end
        if parami == 3
            axis([0 500 0 500])
        end
        plotUnity;
        axis square
    end
end