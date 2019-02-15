function statsGeneral = calcStatsMS_noMotionInfo(tracks,minTrackLen,...
    timeBetweenFrames,removePotArtifacts)
%CALCSTATSMS_NOMOTIONINFO calculates merge/split statistics without distinguishing between different motion types
%
%SYNOPSIS [statsGeneral,motionTypeFrac,probMSIfTypeFull,probMSIfType12,...
%    numEventsPerCatFull,motionRedTypeFrac,probMSIfRedTypeFull,probMSIfRedType12]...
%    = calcStatsMS(tracks,minTrackLen,probDim,diffAnalysisRes,removePotArtifacts)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       timeBetweenFrames: Time between frames (s).
%                    Optional. Devault: 1.
%       removePotArtifacts: 1 to remove potentially artifactual merges and
%                    splits, resulting for instance from detection
%                    artifact, 0 otherwise. 
%                    Optional. Default: 1.
%
%OUTPUT statsGeneral: Row vector with entries: 
%                     (1) Number of features per frame.
%                     (2) Number of tracks with length >= minTrackLen.
%                     (3) Number of tracks segments belonging to the tracks
%                         with length >= minTrackLen.
%                     (4) Probability of a feature merging.
%                     (5) Probability of a feature splitting.
%                     (6) Rate of a feature merging (per unit time).
%                     (7) Rate of a feature splitting (per unit time).
%
%Khuloud Jaqaman, March 2013

%% input

if nargin < 1 || isempty(tracks)
    disp('calcStatsMS: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(timeBetweenFrames)
    timeBetweenFrames = 1;
end

if nargin < 4 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%put tracks in matrix format
[tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);

%get number of track segments
numTrackSegments = size(tracksMat,1);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% merges/splits vs. features

%initialize lists of merges/splits
listOfMerges = [];
listOfSplits = [];

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get the sequence of events of this compound track
    seqOfEvents = tracks(iTrack).seqOfEvents;

    %if requested, remove splits and merges that are most likely artifacts
    if removePotArtifacts
        seqOfEvents = removeSplitMergeArtifacts(seqOfEvents,0);
    end

    %find where this track has merges/splits
    indxMS = find(~isnan(seqOfEvents(:,4)));

    %go over these merges/splits
    for iMS = indxMS'

        %determine whether it's a merge (2) or a split (1)
        msType = seqOfEvents(iMS,2);

        %get the indices of the participating segments within the compound
        %track
        segmentsMS = seqOfEvents(iMS,3:4);

        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;

        %add this merge/split to the list of merges/splits
        if msType == 2 %merge
            listOfMerges = [listOfMerges; segmentsMS]; %#ok<AGROW>
        else %split
            listOfSplits = [listOfSplits; segmentsMS]; %#ok<AGROW>
        end

    end

end

%get total number of merges and splits
numMergesTot = size(listOfMerges,1);
numSplitsTot = size(listOfSplits,1);

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate average number of merges/splits per feature - this is the
%probability of a feature merging/splitting
probFeatMerge = aveMergePerFrame / aveFeatPerFrame;
probFeatSplit = aveSplitPerFrame / aveFeatPerFrame;

%calculate the rate of a feature merging/splitting
rateFeatMerge = probFeatMerge / timeBetweenFrames;
rateFeatSplit = probFeatSplit / timeBetweenFrames;

%% output

%general statistics
statsGeneral = [aveFeatPerFrame numTracks numTrackSegments probFeatMerge ...
    probFeatSplit rateFeatMerge rateFeatSplit];

%% ~~~ the end ~~~

