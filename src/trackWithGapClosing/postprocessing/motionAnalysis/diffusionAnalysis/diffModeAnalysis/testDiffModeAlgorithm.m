function [fracTrueFalsePos0,fracContrRes0,modeParam0,modeParamControl0,...
    fracTrueFalsePos1,fracContrRes1,modeParam1,modeParamControl1] = testDiffModeAlgorithm(...
    modeDiffCoef,fracContrData,diffModeDividerStruct,trackLength,numTracks)

%number of modes
numMode = length(modeDiffCoef);

%generate trajectories
numTracksMode = round(fracContrData*numTracks);
trackCounter = [0 cumsum(numTracksMode)];
tracksFinal0 = repmat(struct('tracksCoordAmpCG',[],'seqOfEvents',[]),trackCounter(end),1);
diffModeSim = [];
for j=1:numMode
    for i=trackCounter(j)+1:trackCounter(j+1)
        traj = brownianMotion(2,modeDiffCoef(j),trackLength,0.1);
        traj = traj(2:10:end,:);
        tmp = zeros(1,trackLength*8);
        tmp(1:8:end) = traj(:,1);
        tmp(2:8:end) = traj(:,2);
        tracksFinal0(i).tracksCoordAmpCG = tmp;
        tracksFinal0(i).seqOfEvents = [1 1 1 NaN; trackLength 2 1 NaN];
    end
    diffModeSim = [diffModeSim; j*ones(numTracksMode(j),1)]; %#ok<AGROW>
end

%from error free data

%get diffusion modes
[modeParam0,~,modeParamControl0] = getDiffModes(tracksFinal0,5,0.01,1,5,2,'test',[],1);

%classify trajectories
diffModeAnalysisRes0 = trackDiffModeAnalysis(tracksFinal0,diffModeDividerStruct);

%analyze classification
diffModeRes0 = vertcat(diffModeAnalysisRes0.diffMode);
fracTrueFalsePos0 = NaN(numMode);
for i=1:numMode
    for j=1:numMode
        fracTrueFalsePos0(j,i) = length(find(diffModeSim==i&diffModeRes0==j))/length(find(diffModeSim==i));
    end
end
fracContrRes0 = hist(diffModeRes0,1:numMode);
fracContrRes0 = fracContrRes0 / sum(fracContrRes0);

%from data with localization error
tracksFinal1 = tracksFinal0;
for i=1:trackCounter(end)
    tracksFinal1(i).tracksCoordAmpCG(1:8:end) = tracksFinal1(i).tracksCoordAmpCG(1:8:end) + randn(1,trackLength)*0.25;
    tracksFinal1(i).tracksCoordAmpCG(2:8:end) = tracksFinal1(i).tracksCoordAmpCG(2:8:end) + randn(1,trackLength)*0.25;
    tracksFinal1(i).tracksCoordAmpCG(5:8:end) = 0.25;
    tracksFinal1(i).tracksCoordAmpCG(6:8:end) = 0.25;
end

%get diffusion modes
[modeParam1,~,modeParamControl1] = getDiffModes(tracksFinal1,5,0.01,1,5,2,'test',[],1);

%classify trajectories
diffModeAnalysisRes1 = trackDiffModeAnalysis(tracksFinal1,diffModeDividerStruct);

%analyze classification
diffModeRes1 = vertcat(diffModeAnalysisRes1.diffMode);
fracTrueFalsePos1 = NaN(numMode);
for i=1:numMode
    for j=1:numMode
        fracTrueFalsePos1(j,i) = length(find(diffModeSim==i&diffModeRes1==j))/length(find(diffModeSim==i));
    end
end
fracContrRes1 = hist(diffModeRes1,1:numMode);
fracContrRes1 = fracContrRes1 / sum(fracContrRes1);
