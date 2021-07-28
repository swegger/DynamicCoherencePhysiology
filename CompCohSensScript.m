%% CompCohSensScript
%
%
%
%%

%% Variables
subject = 'ar';

%% Load dcp objects
load('dcpObjects20210406.mat')

%% Cycle through objects 
unitInd = 0;
for filei = 1:length(dcp)
    disp(['File ' num2str(filei) ' of ' num2str(length(dcp))])
    % Build appropriate initiateCohObj and dynamicCohObj
    dynamicCoh = dynamicCohObj(subject,dcp{filei}.datapath);
    dynamicCoh = assertSpikesExtracted(dynamicCoh,dcp{filei}.spikesExtracted);
    dynamicCoh = unitsIndex(dynamicCoh);
    dynamicCoh = dynamicCohTrials(dynamicCoh,1:5000);
    dynamicCoh = addCoh(dynamicCoh);
    
    initCoh = initiateCohObj(subject,dcp{filei}.datapath);
    initCoh = assertSpikesExtracted(initCoh,dcp{filei}.spikesExtracted);
    initCoh = unitsIndex(initCoh);
    initCoh = initiateCohTrials(initCoh,1:5000);
    
    if isempty(dynamicCoh.sequences)
        disp('No dynamicCoh trials')
    elseif isempty(initCoh.coh)
        disp('No initCoh trials')
    else
        if any(dynamicCoh.unitIndex ~= initCoh.unitIndex)
            disp('Unit numbers disagree')
        else
            for uniti = 1:length(dynamicCoh.unitIndex)
                unitInd = unitInd+1;
                [m(unitInd), ste(unitInd), sensitivity(unitInd)] = ...
                    compareCohTuningInitDynCoh(initCoh,dynamicCoh,...
                    dynamicCoh.unitIndex(uniti));
                
                sens.init(:,:,unitInd) = sensitivity(unitInd).init.theta;
                sens.dyn(:,:,unitInd) = sensitivity(unitInd).dyn.theta;
                
                id(unitInd).file = initCoh.datapath(end-9:end);
                id(unitInd).index = [filei, uniti];
            end
        end
    end
end

%% Compare sensitivities
lim = 0.1;
figure;
subplot(3,3,1)
plot(squeeze(sens.init(1,1,:)),squeeze(sens.init(1,2,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,1,:)),squeeze(sens.init(1,2,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 1')
ylabel('Init 2')


subplot(3,3,2)
plot(squeeze(sens.init(1,1,:)),squeeze(sens.init(1,3,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,1,:)),squeeze(sens.init(1,3,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 1')
ylabel('Init 3')

subplot(3,3,3)
plot(squeeze(sens.init(1,2,:)),squeeze(sens.init(1,3,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,2,:)),squeeze(sens.init(1,3,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 2')
ylabel('Init 3')

subplot(3,3,4)
plot(squeeze(sens.init(1,1,:)),squeeze(sens.dyn(1,2,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,1,:)),squeeze(sens.dyn(1,2,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 1')
ylabel('Dyn 2')

subplot(3,3,5)
plot(squeeze(sens.init(1,2,:)),squeeze(sens.dyn(1,2,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,2,:)),squeeze(sens.dyn(1,2,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 2')
ylabel('Dyn 2')

subplot(3,3,6)
plot(squeeze(sens.dyn(1,3,:)),squeeze(sens.dyn(1,2,:)),'o')
hold on
c = corrcoef(squeeze(sens.dyn(1,3,:)),squeeze(sens.dyn(1,2,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Dyn 3')
ylabel('Dyn 2')

subplot(3,3,7)
plot(squeeze(sens.init(1,1,:)),squeeze(sens.dyn(1,3,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,1,:)),squeeze(sens.dyn(1,3,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 1')
ylabel('Dyn 3')

subplot(3,3,8)
plot(squeeze(sens.init(1,3,:)),squeeze(sens.dyn(1,3,:)),'o')
hold on
c = corrcoef(squeeze(sens.init(1,3,:)),squeeze(sens.dyn(1,3,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Init 3')
ylabel('Dyn 3')

subplot(3,3,9)
plot(squeeze(sens.dyn(1,2,:)),squeeze(sens.dyn(1,3,:)),'o')
hold on
c = corrcoef(squeeze(sens.dyn(1,2,:)),squeeze(sens.dyn(1,3,:)));
text(1.1*lim,0.9*lim,num2str(c(1,2)))
axis([-lim lim -lim lim])
axis square
plotUnity;
plotHorizontal(0);
plotVertical(0);
xlabel('Dyn 2')
ylabel('Dyn 3')

%% Examples

sigLevel = 12;
windows = [50,150; 650,750; 950,1050];
t = -100:1600;
for uniti = 1:length(sensitivity)
    initBIC(uniti,:) = sensitivity(uniti).init.deltaBIC;
    dynBIC(uniti,:) = sensitivity(uniti).dyn.deltaBIC;
    
    initTheta(uniti,:) = sensitivity(uniti).init.theta(1,:);
    dynTheta(uniti,:) = sensitivity(uniti).dyn.theta(1,:);
end

%% 'F' neuron Example
exIndex = [length(dcp)-1, find(dcp{length(dcp)-1}.unitIndex == 140)];
listIndex = find(ismember(vertcat(id(:).index), exIndex, 'rows'));
initCoh = initiateCohObj(subject,dcp{exIndex(1)}.datapath);
initCoh = assertSpikesExtracted(initCoh,dcp{exIndex(1)}.spikesExtracted);
initCoh = unitsIndex(initCoh);
initCoh = initiateCohTrials(initCoh,1:5000);
initCoh = cohConditionedRates(initCoh,'width',50,'dirs',0,'speeds',10);

dynamicCoh = dynamicCohObj(subject,dcp{exIndex(1)}.datapath);
dynamicCoh = assertSpikesExtracted(dynamicCoh,dcp{exIndex(1)}.spikesExtracted);
dynamicCoh = unitsIndex(dynamicCoh);
dynamicCoh = dynamicCohTrials(dynamicCoh,1:5000);
dynamicCoh = dynamicCohSeqConditionedRates(dynamicCoh,'width',50,'dirs',0);
dynamicCoh = addCoh(dynamicCoh);


plotInitSeqComp(initCoh,dynamicCoh,140,windows,sensitivity(listIndex),m(listIndex),ste(listIndex),'t',t)

%% Init1 and Dyn2 Significantly tuned to coherence
sigTunedInit1 = initBIC(:,1) >= 12;
sigTunedDyn2 = dynBIC(:,2) >= 12;

exInds = find(sigTunedInit1 & sigTunedDyn2);

for indi = 1:length(exInds)
    exIndex = id(exInds(indi)).index;
    initCoh = initiateCohObj(subject,dcp{exIndex(1)}.datapath);
    initCoh = assertSpikesExtracted(initCoh,dcp{exIndex(1)}.spikesExtracted);
    initCoh = unitsIndex(initCoh);
    initCoh = initiateCohTrials(initCoh,1:5000);
    initCoh = cohConditionedRates(initCoh,'width',50,'dirs',0,'speeds',10);
    
    dynamicCoh = dynamicCohObj(subject,dcp{exIndex(1)}.datapath);
    dynamicCoh = assertSpikesExtracted(dynamicCoh,dcp{exIndex(1)}.spikesExtracted);
    dynamicCoh = unitsIndex(dynamicCoh);
    dynamicCoh = dynamicCohTrials(dynamicCoh,1:5000);
    dynamicCoh = dynamicCohSeqConditionedRates(dynamicCoh,'width',50,'dirs',0);
    dynamicCoh = addCoh(dynamicCoh);
    
    plotInitSeqComp(initCoh,dynamicCoh,initCoh.unitIndex(exIndex(2)),windows,...
        sensitivity(exInds(indi)),m(exInds(indi)),ste(exInds(indi)),'t',t)
end

%% Init1 significantly tuned, Dyn2 opposite tuning to coherence
exInds = find(initTheta(:,1).*dynTheta(:,2) < 0 & initBIC(:,1) > 10);

for indi = 1:length(exInds)
    exIndex = id(exInds(indi)).index;
    initCoh = initiateCohObj(subject,dcp{exIndex(1)}.datapath);
    initCoh = assertSpikesExtracted(initCoh,dcp{exIndex(1)}.spikesExtracted);
    initCoh = unitsIndex(initCoh);
    initCoh = initiateCohTrials(initCoh,1:5000);
    initCoh = cohConditionedRates(initCoh,'width',50,'dirs',0,'speeds',10);
    
    dynamicCoh = dynamicCohObj(subject,dcp{exIndex(1)}.datapath);
    dynamicCoh = assertSpikesExtracted(dynamicCoh,dcp{exIndex(1)}.spikesExtracted);
    dynamicCoh = unitsIndex(dynamicCoh);
    dynamicCoh = dynamicCohTrials(dynamicCoh,1:5000);
    dynamicCoh = dynamicCohSeqConditionedRates(dynamicCoh,'width',50,'dirs',0);
    dynamicCoh = addCoh(dynamicCoh);
    
    plotInitSeqComp(initCoh,dynamicCoh,initCoh.unitIndex(exIndex(2)),windows,...
        sensitivity(exInds(indi)),m(exInds(indi)),ste(exInds(indi)),'t',t)
end

%% All the figures
saveOpts.On = true;
filenameBase = ['/home/seth/Projects/DynamicCoherencePhysiology/' dcp{1}.sname '/initCohDynCohComp/All/'];

for indi = 1:size(initTheta,1)
    exIndex = id(indi).index;
    initCoh = initiateCohObj(subject,dcp{exIndex(1)}.datapath);
    initCoh = assertSpikesExtracted(initCoh,dcp{exIndex(1)}.spikesExtracted);
    initCoh = unitsIndex(initCoh);
    initCoh = initiateCohTrials(initCoh,1:5000);
    initCoh = cohConditionedRates(initCoh,'width',50,'dirs',0,'speeds',10);
    
    dynamicCoh = dynamicCohObj(subject,dcp{exIndex(1)}.datapath);
    dynamicCoh = assertSpikesExtracted(dynamicCoh,dcp{exIndex(1)}.spikesExtracted);
    dynamicCoh = unitsIndex(dynamicCoh);
    dynamicCoh = dynamicCohTrials(dynamicCoh,1:5000);
    dynamicCoh = dynamicCohSeqConditionedRates(dynamicCoh,'width',50,'dirs',0);
    dynamicCoh = addCoh(dynamicCoh);
    
    saveOpts.filename = [filenameBase dcp{1}.sname initCoh.datapath(end-6:end) '_' num2str(initCoh.unitIndex(exIndex(2)))];
    plotInitSeqComp(initCoh,dynamicCoh,initCoh.unitIndex(exIndex(2)),windows,...
        sensitivity(indi),m(indi),ste(indi),'t',t,...
        'closeFig',true,'saveOpts',saveOpts)
end

