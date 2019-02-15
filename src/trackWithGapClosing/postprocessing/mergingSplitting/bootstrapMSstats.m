function msBootstrapRes = bootstrapMSstats(tracksAll,minTrackLen,probDim,...
    diffAnalysisAll,numRep,removePotArtifacts)

%Output: msBootstrapRes contains 16 columns, storing;
%        (1) overall probability of merging/splitting.
%        (2-5) probability of undergoing linear, isotropic-unconfined,
%        isotropic-confined, isotropic-undetermine and completely
%        undetermined motion.
%        (6-11) probability of merging/splitting if undergoing the five
%        motion types listed above, as calculated by method 1 (more dynamic
%        motion type prevails).
%        (12-16) probability of merging/splitting if undergoing the five
%        motino types lister above, as calculated by method 2 (no
%        classification of more or less dynamic).


%remove tracks shorter than minTrackLen
criteria.lifeTime.min = minTrackLen;
indxKeep = chooseTracks(tracksAll,criteria);
tracksAll = tracksAll(indxKeep);
diffAnalysisAll = diffAnalysisAll(indxKeep);

%get number of available tracks
numTracks = length(tracksAll);

%reserve memory for results
msBootstrapRes = NaN(numRep,16);

%repeat numRep times
for iRep = 1 : numRep

    %select a random subset, with replacement
    if iRep ~= 1
        randIndx = randsample((1:numTracks),numTracks,'true');
    else
        randIndx = 1:numTracks;
    end
    tracksSample = tracksAll(randIndx);
    diffAnalysisSample = diffAnalysisAll(randIndx);
    
    %analyze MS stats
    [statsGeneral,motionTypeFrac,~,probMSIfType12] = calcStatsMS(tracksSample,...
        minTrackLen,probDim,diffAnalysisSample,removePotArtifacts);
    
    %extract results
    probMS = mean(statsGeneral(4:5)); %overall merging and splitting probabilities
    probType = motionTypeFrac(:,2)'; %overall type probabilities
    prob1 = probMSIfType12(:,1,3)'; %conditional probabilities as calculated in approach 1
    prob2 = probMSIfType12(:,2,3)'; %conditional probabilities as calculated in approach 2

    %store results
    msBootstrapRes(iRep,:) = [probMS probType prob1 prob2];

end

