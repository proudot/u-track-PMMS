function [statsGeneral,motionTypeFrac,probMSIfTypeFull,probMSIfType12,...
    numEventsPerCatFull,motionRedTypeFrac,probMSIfRedTypeFull,probMSIfRedType12]...
    = calcStatsMS(tracks,minTrackLen,probDim,diffAnalysisRes,removePotArtifacts)
%CALCSTATSMS calculates merge/split statistics for linear, Brownian and confined Brownian tracks
%
%SYNOPSIS [statsGeneral,motionTypeFrac,probMSIfTypeFull,probMSIfType12,...
%    numEventsPerCatFull,motionRedTypeFrac,probMSIfRedTypeFull,probMSIfRedType12]...
%    = calcStatsMS(tracks,minTrackLen,probDim,diffAnalysisRes,removePotArtifacts)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffAnalysis1). Optional. If not input, it will
%                    be calculated.
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
%    All the following variables have 5 rows, corresponding to the motion
%    types:
%                     (1) linear, (2) Brownian, (3) confined Brownian,
%                     (4) undetermined >= 5 frames (i.e. not linear but
%                     undetermine whether Brownian or confined), and
%                     (5) undetermined < 5 frames (i.e. can be anything).
%       motionTypeFrac: 5-by-2 array. Columns corrrespond to:
%                     (1) Fraction of track segments in each motion type.
%                     (2) Probability of a feature undergoing a motion type.
%       probMSIfTypeFull: 5-by-5-by-3 array showing double conditional
%                     probability of merging/splitting for the different
%                     motion type combinations. Columns are transpose of
%                     rows. Third dimension corresponds to (1) merges, (2)
%                     splits and (3) their average.
%       probMSIfType12: 5-by-2-by-3 array. Columns correspond to
%                     (1) Conditional probability of a feature
%                         merging/splitting IF undergoing a motion type,
%                         where the more dynamic motion type is dominant.
%                     (2) Conditional probability of a feature
%                         merging/splitting IF undergoing a motion type,
%                         regardless of the other feature's motion type.
%                     Third dimension corresponds to (1) merges, (2) splits
%                     and (3) their average.
%       numEventsPerCatFull: 5-by-5-by-3 array showing number of
%                     merging/splitting events observed per motion type.
%                     Rows and columns correspond to motion types. Third
%                     dimension corresponds to (1) merging, (2) splitting
%                     and (3) sum of merging and splitting events.
%    All the following variables have 3 rows, corresponding to the motion
%    types:
%                     (1) linear, (2) not linear, i.e. sum of 2-4 above,
%                     and (3) undetermined < 5 frames.
%       motionRedTypeFrac: 3-by-2 array. Columns corrrespond to:
%                     (1) Fraction of track segments in each motion type.
%                     (2) Probability of a feature undergoing a motion type.
%       probMSIfRedTypeFull: 3-by-3-by-3 array showing double conditional
%                     probability of merging/splitting for the different
%                     motion type combinations. Columns are transpose of
%                     rows. Third dimension corresponds to (1) merges, (2)
%                     splits and (3) their average.
%       probMSIfRedType12: 3-by-2-by-3 array. Columns correspond to
%                     (1) Conditional probability of a feature
%                         merging/splitting IF undergoing a motion type,
%                         where the more dynamic motion type is dominant.
%                     (2) Conditional probability of a feature
%                         merging/splitting IF undergoing a motion type,
%                         regardless of the other feature's motion type.
%                     Third dimension corresponds to (1) merges, (2) splits
%                     and (3) their average.
%
%Khuloud Jaqaman, December 2007

%% input

if nargin < 1 || isempty(tracks)
    disp('calcStatsMS: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(diffAnalysisRes)
    [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,1,probDim,...
        1,[0.05 0.1],0);
    if errFlag
        return
    end
end

if nargin < 5 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);
diffAnalysisRes = diffAnalysisRes(indx);

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

%% track segment types

%get track segment classification from diffusion analysis
trackSegmentClass = vertcat(diffAnalysisRes.classification);

%get track segment length
trackSegmentLength = getTrackSEL(tracksMat);
trackSegmentLength = trackSegmentLength(:,3);

%get indices of linear, Brownian, confined Brownian and undetermined tracks
indxLin    = find( trackSegmentClass(:,1) == 1 | trackSegmentClass(:,2) == 3 );
indxBrown  = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 2 );
indxConf   = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 1 );
indxUndet1 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
    & trackSegmentLength >= 5);
indxUndet2 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
    & trackSegmentLength < 5);

%store track segment type in an array
%-1 = undetermined & track length < 5, 0 = undetermined & track length >= 5, 
%1 = confined Brownian, 2 = Brownian, 3 = linear
trackSegmentType = zeros(numTrackSegments,1);
trackSegmentType(indxUndet2) = -1;
trackSegmentType(indxConf)   = 1;
trackSegmentType(indxBrown)  = 2;
trackSegmentType(indxLin)    = 3;

%calculate number of track segments per category
numSegmentsType = [length(indxLin) length(indxBrown) length(indxConf) ...
    length(indxUndet1) length(indxUndet2)]';

%calculate number of track segments per reduced category (i.e. 0, 1 and 2
%are lumped together into one non-linear category)
numSegmentsRedType = [numSegmentsType(1) sum(numSegmentsType(2:4)) numSegmentsType(5)]';

%calculate fraction of track segments falling in each category
fracSegmentsType = numSegmentsType / numTrackSegments;

%calculate fraction of track segments falling in each reduced category
fracSegmentsRedType = numSegmentsRedType / numTrackSegments;

%calculate number of features in each category and each reduced category
numFeatType = [length(find(tracksIndxMat(indxLin,:))) ...
    length(find(tracksIndxMat(indxBrown,:))) ...
    length(find(tracksIndxMat(indxConf,:))) ...
    length(find(tracksIndxMat(indxUndet1,:))) ...
    length(find(tracksIndxMat(indxUndet2,:)))]';
numFeatRedType = [numFeatType(1) sum(numFeatType(2:4)) numFeatType(5)]';

%get fraction of features in each category - this is the probability of a
%feature undergoing a certain motion type
probFeatType = numFeatType / numFeatTot;

%repeat for reduced categories
probFeatRedType = numFeatRedType / numFeatTot;

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

%% merges/splits vs. type

%get the types of segments participating in merges and splits
if numMergesTot > 0
    listOfMergeTypes = [trackSegmentType(listOfMerges(:,1)) ...
        trackSegmentType(listOfMerges(:,2))];
else
    listOfMergeTypes = zeros(0,2);
end
if numSplitsTot > 0
    listOfSplitTypes = [trackSegmentType(listOfSplits(:,1)) ...
        trackSegmentType(listOfSplits(:,2))];
else
    listOfSplitTypes = zeros(0,2);
end

%sort the lists so that the "larger" (more dynamic) type is in the first
%column
listOfMergeTypes = sort(listOfMergeTypes,2,'descend');
listOfSplitTypes = sort(listOfSplitTypes,2,'descend');

%get number of merges/splits based on the types of the two track segments
%involved
typeList = [3 2 1 0 -1];
[numMergesTypeType1,numMergesTypeType2] = deal(NaN(5));
for iType = 1 : 5
    for jType = iType : 5
        [numMergesTypeType1(iType,jType),numMergesTypeType2(iType,jType),...
            numMergesTypeType2(jType,iType)] = deal(length(find( ...
            listOfMergeTypes(:,1)==typeList(iType) & ...
            listOfMergeTypes(:,2)==typeList(jType) )));
    end
end
[numSplitsTypeType1,numSplitsTypeType2] = deal(NaN(5));
for iType = 1 : 5
    for jType = iType : 5
        [numSplitsTypeType1(iType,jType),numSplitsTypeType2(iType,jType),...
            numSplitsTypeType2(jType,iType)] = deal(length(find( ...
            listOfSplitTypes(:,1)==typeList(iType) & ...
            listOfSplitTypes(:,2)==typeList(jType) )));
    end
end

%simplification 1: 
%classify a merge/split type based on the more dynamic of its two segments
%(linear > Brownian > confined Browian > (undetermined&length>=5) > (undetermined&length<5))
%for example, if a linear track segment merges with a Brownian track
%segment, the merge is classified as linear
numMergesType1 = nansum(numMergesTypeType1,2);
numSplitsType1 = nansum(numSplitsTypeType1,2);

%simplification 2:
%classify a merge/split type based on one of its two segments regardless of
%the other segment, i.e. no ranking of dynamicity in this case
numMergesType2 = sum(numMergesTypeType2,2);
numSplitsType2 = sum(numSplitsTypeType2,2);

%also get number of merges/splits based on reduced track segment types
%i.e. 2, 1 and 0 are lumped into one non-linear type
numMergesRedTypeType1 = NaN(3);
numMergesRedTypeType1(1,:) = [numMergesTypeType1(1,1) sum(numMergesTypeType1(1,2:4)) numMergesTypeType1(1,5)];
numMergesRedTypeType1(2,2:3) = [nansum(nansum(numMergesTypeType1(2:4,2:4))) sum(numMergesTypeType1(2:4,5))];
numMergesRedTypeType1(3,3) = numMergesTypeType1(5,5);
tmpMat(:,:,1) = numMergesRedTypeType1;
tmpMat(:,:,2) = numMergesRedTypeType1';
numMergesRedTypeType2 = max(tmpMat,[],3);
numSplitsRedTypeType1 = NaN(3);
numSplitsRedTypeType1(1,:) = [numSplitsTypeType1(1,1) sum(numSplitsTypeType1(1,2:4)) numSplitsTypeType1(1,5)];
numSplitsRedTypeType1(2,2:3) = [nansum(nansum(numSplitsTypeType1(2:4,2:4))) sum(numSplitsTypeType1(2:4,5))];
numSplitsRedTypeType1(3,3) = numSplitsTypeType1(5,5);
tmpMat(:,:,1) = numSplitsRedTypeType1;
tmpMat(:,:,2) = numSplitsRedTypeType1';
numSplitsRedTypeType2 = max(tmpMat,[],3);

%simplification 1:
numMergesRedType1 = nansum(numMergesRedTypeType1,2);
numSplitsRedType1 = nansum(numSplitsRedTypeType1,2);

%simplification 2:
numMergesRedType2 = sum(numMergesRedTypeType2,2);
numSplitsRedType2 = sum(numSplitsRedTypeType2,2);

%get the fraction of merge/split types - this is the conditional
%probability of having a certain motion type if merging/splitting

%simplification 1:
probTypeIfMerge1 = numMergesType1 / numMergesTot;
probTypeIfSplit1 = numSplitsType1 / numSplitsTot;

%simplification 2:
probTypeIfMerge2 = numMergesType2 / numMergesTot;
probTypeIfSplit2 = numSplitsType2 / numSplitsTot;

%full:
probTypeTypeIfMerge = numMergesTypeType1 / numMergesTot;
probTypeTypeIfSplit = numSplitsTypeType1 / numSplitsTot;

%reduced types, simplification 1:
probRedTypeIfMerge1 = numMergesRedType1 / numMergesTot;
probRedTypeIfSplit1 = numSplitsRedType1 / numSplitsTot;

%reduced types, simplification 2:
probRedTypeIfMerge2 = numMergesRedType2 / numMergesTot;
probRedTypeIfSplit2 = numSplitsRedType2 / numSplitsTot;

%reduced types, full:
probRedTypeTypeIfMerge = numMergesRedTypeType1 / numMergesTot;
probRedTypeTypeIfSplit = numSplitsRedTypeType1 / numSplitsTot;

%calculate the conditional probability of a feature undergoing a
%merge/split IF undergoing a certain motion type

%simplification 1:
probMergeIfType1 = probTypeIfMerge1 * probFeatMerge ./ probFeatType;
probSplitIfType1 = probTypeIfSplit1 * probFeatSplit ./ probFeatType;

%simplification 2:
probMergeIfType2 = probTypeIfMerge2 * probFeatMerge ./ probFeatType;
probSplitIfType2 = probTypeIfSplit2 * probFeatSplit ./ probFeatType;

%full:
probFeatTypeType = probFeatType * probFeatType';
probMergeIfTypeType = probTypeTypeIfMerge * probFeatMerge ./ probFeatTypeType;
probSplitIfTypeType = probTypeTypeIfSplit * probFeatSplit ./ probFeatTypeType;

%reduced types, simplification 1:
probMergeIfRedType1 = probRedTypeIfMerge1 * probFeatMerge ./ probFeatRedType;
probSplitIfRedType1 = probRedTypeIfSplit1 * probFeatSplit ./ probFeatRedType;

%reduced types, simplification 2:
probMergeIfRedType2 = probRedTypeIfMerge2 * probFeatMerge ./ probFeatRedType;
probSplitIfRedType2 = probRedTypeIfSplit2 * probFeatSplit ./ probFeatRedType;

%reduced types, full:
probFeatRedTypeType = probFeatRedType * probFeatRedType';
probMergeIfRedTypeType = probRedTypeTypeIfMerge * probFeatMerge ./ probFeatRedTypeType;
probSplitIfRedTypeType = probRedTypeTypeIfSplit * probFeatSplit ./ probFeatRedTypeType;

%% output

%general statistics
statsGeneral = [aveFeatPerFrame numTracks numTrackSegments probFeatMerge probFeatSplit];

%motion type probabilities
%full and reduced
motionTypeFrac = [fracSegmentsType probFeatType];
motionRedTypeFrac = [fracSegmentsRedType probFeatRedType];

%merging and splitting probability per motion type
%full
probMSIfTypeFull(:,:,1) = probMergeIfTypeType;
probMSIfTypeFull(:,:,2) = probSplitIfTypeType;
probMSIfTypeFull(:,:,3) = nanmean(probMSIfTypeFull(:,:,1:2),3);
%simplified
probMSIfType12(:,:,1) = [probMergeIfType1 probMergeIfType2];
probMSIfType12(:,:,2) = [probSplitIfType1 probSplitIfType2];
probMSIfType12(:,:,3) = nanmean(probMSIfType12(:,:,1:2),3);

%merging and splitting probability per reduced motion type
%full
probMSIfRedTypeFull(:,:,1) = probMergeIfRedTypeType;
probMSIfRedTypeFull(:,:,2) = probSplitIfRedTypeType;
probMSIfRedTypeFull(:,:,3) = nanmean(probMSIfRedTypeFull(:,:,1:2),3);
%simplified
probMSIfRedType12(:,:,1) = [probMergeIfRedType1 probMergeIfRedType2];
probMSIfRedType12(:,:,2) = [probSplitIfRedType1 probSplitIfRedType2];
probMSIfRedType12(:,:,3) = nanmean(probMSIfRedType12(:,:,1:2),3);

%number of events per category
%full
numEventsPerCatFull(:,:,1) = numMergesTypeType1;
numEventsPerCatFull(:,:,2) = numSplitsTypeType1;
numEventsPerCatFull(:,:,3) = sum(numEventsPerCatFull(:,:,1:2),3);

%% ~~~ the end ~~~

