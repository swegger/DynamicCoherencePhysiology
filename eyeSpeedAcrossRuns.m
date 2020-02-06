function [mEye, stdEye, nEye] = eyeSpeedAcrossRuns(dynamicCoh)
%%
%
%
%%

dirs = [];
speeds = [];
sequences = [];
eyeSpeed = [];
for i = 1:length(dynamicCoh)
    
    if ~isempty(dynamicCoh{i}.directions)
        dirs = [dirs; dynamicCoh{i}.directions];
        speeds = [speeds; dynamicCoh{i}.directions];
        sequences = [sequences; dynamicCoh{i}.sequences];
        for ti = 1:length(dynamicCoh{i}.eye)
            eyeSpeed = [eyeSpeed; sqrt(dynamicCoh{i}.eye(ti).hvel.^2 + dynamicCoh{i}.eye(ti).vvel.^2)];
        end
    end
end

%%
seqs = unique(sequences);
dir = 0;
for seqi = 1:length(seqs)
    dirvec = dirs == dir;
    seqvec = sequences == seqs(seqi);
    
    mEye(:,seqi) = nanmean(eyeSpeed(seqvec & dirvec,:),1);
    stdEye(:,seqi) = nanstd(eyeSpeed(seqvec & dirvec,:),1);
    nEye(:,seqi) = sum(seqvec & dirvec);
end

%% Plot
colors = colormap('lines');
t = 0:size(mEye,1)-1;

figure('Position',[752 580 804 525])
subplot(2,1,1)
for seqi = 1:length(seqs)    
    mEd = mEye(:,seqi);
    steEd = stdEye(:,seqi)./sqrt(nEye(:,seqi));
    
    
    patchProps.FaceAlpha = 0.3;
    patchProps.FaceColor = colors(seqi,:);
    myPatch(t(:),mEd(:),steEd(:),'patchProperties',patchProps);
    hold on
end
for seqi = 1:length(seqs)
    plot(t,mEye(:,seqi),'Color',colors(seqi,:))
    hold on
end
plotVertical([150 150+0:300:1500]);
xlabel('Time (ms)')
ylabel('Mean eye speed')
mymakeaxis(gca)

controlSeq = 5;
subplot(2,1,2)
for seqi = 1:length(seqs)    
    mEd = mEye(:,seqi)-mEye(:,controlSeq);
    steEd = sqrt((stdEye(:,seqi).^2+stdEye(:,seqi).^2)./(nEye(:,seqi)+nEye(:,controlSeq)));
    
    
    patchProps.FaceAlpha = 0.3;
    patchProps.FaceColor = colors(seqi,:);
    myPatch(t(:),mEd(:),steEd(:),'patchProperties',patchProps);
    hold on
end
for seqi = 1:length(seqs)
    plot(t,mEye(:,seqi)-mEye(:,controlSeq))
    hold on
end
plotVertical([150 150+0:300:1500]);
xlabel('Time (ms)')
ylabel('Mean eye minus control')
mymakeaxis(gca)