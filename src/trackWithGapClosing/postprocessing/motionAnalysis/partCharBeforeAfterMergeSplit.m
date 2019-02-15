function [charBeforeAfterMerge,charBeforeAfterSplit] = ...
    partCharBeforeAfterMergeSplit(tracks,minTrackLen,probDim,...
    diffAnalysisRes,removePotArtifacts)
%partCharBeforeAfterMergeSplit calculates particle properties before and after merges and splits
%
%SYNOPSIS [charBeforeAfterMerge,charBeforeAfterSplit] = ...
%    partCharBeforeAfterMergeSplit(tracks,minTrackLen,probDim,...
%    diffAnalysisRes,removePotArtifacts)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used.
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
%OUTPUT ...
%
%Khuloud Jaqaman, October 2010

%% Input

if nargin < 1 || isempty(tracks)
    disp('partCharBeforeAfterMergeSplit: Missing input argument!');
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
        1,[0.05 0.05],0);
    if errFlag
        return
    end
end

if nargin < 5 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

%% Preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);
diffAnalysisRes = diffAnalysisRes(indx);

%get number of tracks
numTracks = length(tracks);

%put tracks in matrix format
[tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);
xCoordMat = tracksMat(:,1:8:end);
yCoordMat = tracksMat(:,2:8:end);
zCoordMat = tracksMat(:,3:8:end);
ampMat = tracksMat(:,4:8:end);

%get number of track segments and number of frames
[numTrackSegments,numFrames] = size(ampMat);

%% Track segment types

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

%% Merge and split characteristics

%initialize variable for storage
infoBeforeAfter = [];

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
        
        %get the merge/split time
        msTime = seqOfEvents(iMS,1);
        
        %determine whether it's a merge (2) or a split (1)
        msType = seqOfEvents(iMS,2);
        
        %get the indices of the participating segments within the compound
        %track
        segmentsMS = seqOfEvents(iMS,3:4);

        %get their indices in the global segment matrix
        segmentsMS = segmentsMS + trackStartRow(iTrack) - 1;
        
        %get the motion type
        motionType = max(trackSegmentType(segmentsMS));
        
        %calculate some indices taking into account movie start and end
        %times
        msTimeMinus10 = max(msTime-10,1);
        msTimeMinus1 = max(msTime-1,1);
        msTimePlus9 = min(msTime+9,numFrames);
        
        %calculate the intensity characteristics before and after
        %the merge or the split
        intBefore = ampMat(segmentsMS,msTimeMinus10:msTimeMinus1);
        intBefore = nanmean(intBefore(:));
        intAfter  = ampMat(segmentsMS,msTime:msTimePlus9);
        intAfter  = nanmean(intAfter(:));
        
        %calculate the displacement characteristics before the merge or the
        %split
        xBefore = xCoordMat(segmentsMS,msTimeMinus10:msTimeMinus1);
        yBefore = yCoordMat(segmentsMS,msTimeMinus10:msTimeMinus1);
        zBefore = zCoordMat(segmentsMS,msTimeMinus10:msTimeMinus1);
        dispBefore = sqrt(diff(xBefore,[],2).^2 + diff(yBefore,[],2).^2 + ...
            diff(zBefore,[],2).^2);
        dispBefore = nanmean(dispBefore(:));

        %calculate the displacement characteristics after the merge or the
        %split
        xAfter  = xCoordMat(segmentsMS,msTime:msTimePlus9);
        yAfter  = yCoordMat(segmentsMS,msTime:msTimePlus9);
        zAfter  = zCoordMat(segmentsMS,msTime:msTimePlus9);
        dispAfter = sqrt(diff(xAfter,[],2).^2 + diff(yAfter,[],2).^2 + ...
            diff(zAfter,[],2).^2);
        dispAfter = nanmean(dispAfter(:));
        
        %save this information
        infoBeforeAfter = [infoBeforeAfter; [msType motionType intBefore intAfter dispBefore dispAfter]]; %#ok<AGROW>
        
    end

end

%separate merges from splits
charBeforeAfterMerge = infoBeforeAfter(infoBeforeAfter(:,1)==2,2:end);
charBeforeAfterSplit = infoBeforeAfter(infoBeforeAfter(:,1)==1,2:end);


%% ~~~ the end ~~~

