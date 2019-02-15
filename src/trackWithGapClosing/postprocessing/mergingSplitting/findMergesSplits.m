function [mergesInfo,splitsInfo,mergesInfoSpace,splitsInfoSpace,...
    mergesInfoSpaceExt,splitsInfoSpaceExt,mergesInfoAmp,splitsInfoAmp] = ...
    findMergesSplits(tracks,probDim,removePotArtifacts,plotRes,calcTrackType)
%FINDMERGESSPLITS finds the merges and splits in each track and gives back their time and location
%
%SYNOPSIS [mergesInfo,splitsInfo,mergesInfoSpace,splitsInfoSpace,...
%    mergesInfoSpaceExt,splitsInfoSpaceExt] = findMergesSplits(tracks,...
%    probDim,removePotArtifacts,plotRes,calcTrackType)
%
%INPUT  tracks    : Output of trackCloseGapsKalman:
%                   Structure array with number of entries equal to
%                   the number of tracks (or compound tracks when
%                   merging/splitting are considered). Contains the
%                   fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       probDim   : 2 for 2D, 3 for 3D. Optional. Default: 2.
%       removePotArtifacts: 1 to remove potentially artifactual merges and
%                   splits, resulting for instance from detection
%                   artifacts, 0 otherwise. 
%                   Optional. Default: 1.
%       plotRes   : 0 to not plot anything, 1 to make a spatial map of
%                   merges and splits.
%                   Optional. Default: 0.
%       calcTrackType: Estimate track type based on asymmetry. 1 to
%                   estimate, 0 to not estimate. 
%                   Optional. Default: 1.
%
%OUTPUT mergesInfo     : 2D array where first column indicates track number,
%                        second column indicates track type (1 linear, 0 o.w.),
%                        third column indicates number of merges, and 
%                        subsequent columns indicate merge times.
%                        Track type will be 0 if calcTrackType = 0.
%       splitsInfo     : 2D array where first column indicates track number,
%                        second column indicates track type (1 linear, 0 o.w.),
%                        third column indicates number of splits, and 
%                        subsequence columns indicate split times.
%                        Track type will be 0 if calcTrackType = 0.
%       mergesInfoSpace: 2D array that is a continuation of mergesInfoTime,
%                        storing the (x,y,[z])-coordinates of each merge.
%                        Every row corresponds to the same row in
%                        mergesInfo. Every merge gets 2 (in 2D) or 3 (in
%                        3D) columns for x, y and z (if 3D).
%       splitsInfoSpace: 2D array that is a continuation of splitsInfoTime,
%                        storing the (x,y,[z])-coordinates of each split.
%                        Every row corresponds to the same row in
%                        splitsInfo. Every split gets 2 (in 2D) or 3 (in
%                        3D) columns for x, y and z (if 3D).
%       mergesInfoSpaceExt: 3D array that is a continuation of mergesInfoSpace
%                        but replicated in the 3rd dimension, storing the
%                        (x,y,[z])-coordinates of the two particles before
%                        merging. First index of 3rd dimension is for
%                        continuing track segment, second index is for
%                        terminatating track segment.
%       splitsInfoSpaceExt: 3D array that is a continuation of splitsInfoSpace
%                        but replicated in the 3rd dimension, storing the
%                        (x,y,[z])-coordinates of the two particles after
%                        merging. First index of 3rd dimension is for
%                        continuing track segment, second index is for
%                        initiating track segment.
%       mergesInfoAmp  : 2D array that is a continuation of mergesInfo,
%                        with amplitude information stored as amplitude
%                        after merging, amplitude of continuing track 
%                        segment before merging, and amplitude of
%                        terminated track segment before merging.
%       splitsInfoAmp  : 2D array that is a continuation of splitsInfo,
%                        with amplitude information stored as amplitude
%                        before splitting, amplitude of continuing
%                        track segment after splitting, and amplitude of 
%                        starting track segment after splitting.
%
%REMARKS Plotting implemented for 2D only.
%
%Khuloud Jaqaman, October 2007, January 2015

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--findMergesSplits: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

if nargin < 4 || isempty(plotRes)
    plotRes = 0;
end
if probDim ~= 2
    plotRes = 0;
end

if nargin < 5 || isempty(calcTrackType)
    calcTrackType = 1;
end

%get number of tracks
numTracks = length(tracks);

%get number of segments per track
numSegments = getNumSegmentsPerTrack(tracks);

%estimate track types
if calcTrackType
    trackType = getTrackType(tracks,probDim);
else
    trackType = zeros(numTracks,1);
end

%% Merge/split statistics

[mergesInfo,splitsInfo] = deal(zeros(numTracks,max(numSegments)));
[mergesInfoSpace,splitsInfoSpace,mergesInfoAmp,splitsInfoAmp] = deal(repmat(mergesInfo,1,3));
[mergesInfoSpaceExt,splitsInfoSpaceExt] = deal(repmat(mergesInfoSpace,[1 1 2]));

%go over all tracks ...
for iTrack = 1 : numTracks
    
    %get track's sequence of events
    seqOfEvents = tracks(iTrack).seqOfEvents;
    
    %remove splits and merges that are most likely artifacts
    if removePotArtifacts
        seqOfEvents = removeSplitMergeArtifacts(seqOfEvents,0);
    end

    %get track's coordinates and amplitudes
    trackCoordX = tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
    trackCoordY = tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
    trackCoordZ = tracks(iTrack).tracksCoordAmpCG(:,3:8:end);
    trackAmp = tracks(iTrack).tracksCoordAmpCG(:,4:8:end);
    
    %% merges
    
    %find rows with merging information
    indxMerge = find( seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)) );
    
    %get the merge times (absolute and relative to track start)
    mergeTimes = seqOfEvents(indxMerge,1);
    relMergeTimes = mergeTimes - seqOfEvents(1,1) + 1;
    
    %get the track segments that are merged with
    mergeSegment = seqOfEvents(indxMerge,4);
    
    %get the track segments that terminate by merging
    terminateSegment = seqOfEvents(indxMerge,3);
    
    %get the coordinates of each merge
    mergeCoords = [];
    for iMerge = 1 : length(indxMerge)
        mergeCoords = [mergeCoords ...
            trackCoordX(mergeSegment(iMerge),relMergeTimes(iMerge)) ...
            trackCoordY(mergeSegment(iMerge),relMergeTimes(iMerge)) ...
            trackCoordZ(mergeSegment(iMerge),relMergeTimes(iMerge))]; %#ok<AGROW>
    end
    
    %get the coordinates before merging
    mergeCoords21 = [];
    mergeCoords22 = [];
    for iMerge = 1 : length(indxMerge)
        mergeCoords21 = [mergeCoords21 ...
            trackCoordX(mergeSegment(iMerge),relMergeTimes(iMerge)-1) ...
            trackCoordY(mergeSegment(iMerge),relMergeTimes(iMerge)-1) ...
            trackCoordZ(mergeSegment(iMerge),relMergeTimes(iMerge)-1)]; %#ok<AGROW>
        mergeCoords22 = [mergeCoords22 ...
            trackCoordX(terminateSegment(iMerge),relMergeTimes(iMerge)-1) ...
            trackCoordY(terminateSegment(iMerge),relMergeTimes(iMerge)-1) ...
            trackCoordZ(terminateSegment(iMerge),relMergeTimes(iMerge)-1)]; %#ok<AGROW>
    end
    
    %get the amplitudes before and after
    mergeAmps = [];
    for iMerge = 1 : length(indxMerge)
        mergeAmps = [mergeAmps ...
            trackAmp(mergeSegment(iMerge),relMergeTimes(iMerge)) ...
            trackAmp(mergeSegment(iMerge),relMergeTimes(iMerge)-1) ...
            trackAmp(terminateSegment(iMerge),relMergeTimes(iMerge)-1)]; %#ok<AGROW>
    end
    
    %store the merge information for this track
    mergesInfo(iTrack,1:length(mergeTimes)+2) = [trackType(iTrack) ...
        length(mergeTimes) mergeTimes'];
    if ~isempty(mergeTimes)
        mergesInfoSpace(iTrack,1:3*length(mergeTimes)) = mergeCoords;
        mergesInfoSpaceExt(iTrack,1:3*length(mergeTimes),1) = mergeCoords21;
        mergesInfoSpaceExt(iTrack,1:3*length(mergeTimes),2) = mergeCoords22;
        mergesInfoAmp(iTrack,1:3*length(mergeTimes)) = mergeAmps;
    end
        
    %% splits
    
    %find rows with splitting information
    indxSplit = find( seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)) );
    
    %get split times (absolute and relative to track start)
    splitTimes = seqOfEvents(indxSplit,1);
    relSplitTimes = splitTimes - seqOfEvents(1,1) + 1;

    %get the track segments that are split from
    splitSegment = seqOfEvents(indxSplit,4);
    
    %get the track segments that start by splitting
    startSegment = seqOfEvents(indxSplit,3);
    
    %get the coordinates of each split
    splitCoords = [];
    for iSplit = 1 : length(indxSplit)
        splitCoords = [splitCoords ...
            trackCoordX(splitSegment(iSplit),relSplitTimes(iSplit)-1) ...
            trackCoordY(splitSegment(iSplit),relSplitTimes(iSplit)-1) ...
            trackCoordZ(splitSegment(iSplit),relSplitTimes(iSplit)-1)]; %#ok<AGROW>
    end
    
    %get the coordinates after splitting
    splitCoords21 = [];
    splitCoords22 = [];
    for iSplit = 1 : length(indxSplit)
        splitCoords21 = [splitCoords21 ...
            trackCoordX(splitSegment(iSplit),relSplitTimes(iSplit)) ...
            trackCoordY(splitSegment(iSplit),relSplitTimes(iSplit)) ...
            trackCoordZ(splitSegment(iSplit),relSplitTimes(iSplit))]; %#ok<AGROW>
        splitCoords22 = [splitCoords22 ...
            trackCoordX(startSegment(iSplit),relSplitTimes(iSplit)) ...
            trackCoordY(startSegment(iSplit),relSplitTimes(iSplit)) ...
            trackCoordZ(startSegment(iSplit),relSplitTimes(iSplit))]; %#ok<AGROW>
    end
    
    %get the amplitudes before and after
    splitAmps = [];
    for iSplit = 1 : length(indxSplit)
        splitAmps = [splitAmps ...
            trackAmp(splitSegment(iSplit),relSplitTimes(iSplit)-1) ...
            trackAmp(splitSegment(iSplit),relSplitTimes(iSplit)) ...
            trackAmp(startSegment(iSplit),relSplitTimes(iSplit))]; %#ok<AGROW>
    end
    
    %store the split information for this track
    splitsInfo(iTrack,1:length(splitTimes)+2) = [trackType(iTrack) ...
        length(splitTimes) splitTimes'];
    if ~isempty(splitTimes)
        splitsInfoSpace(iTrack,1:3*length(splitTimes)) = splitCoords;
        splitsInfoSpaceExt(iTrack,1:3*length(splitTimes),1) = splitCoords21;
        splitsInfoSpaceExt(iTrack,1:3*length(splitTimes),2) = splitCoords22;
        splitsInfoAmp(iTrack,1:3*length(splitTimes)) = splitAmps;
    end

end

%remove empty columns
if ~isempty(mergesInfo)
    fullIndx = find( sum(mergesInfo(:,2:end))~=0 & ~isnan(sum(mergesInfo(:,2:end))) );
    mergesInfo = [mergesInfo(:,1) mergesInfo(:,1+fullIndx)];
    fullIndx = sum(mergesInfoSpace)~=0 & ~isnan(sum(mergesInfoSpace));
    mergesInfoSpace = mergesInfoSpace(:,fullIndx);
    mergesInfoSpaceExt = mergesInfoSpaceExt(:,fullIndx,:);
    fullIndx = sum(mergesInfoAmp)~=0 & ~isnan(sum(mergesInfoAmp));
    mergesInfoAmp = mergesInfoAmp(:,fullIndx);
end
if ~isempty(splitsInfo)
    fullIndx = find( sum(splitsInfo(:,2:end))~=0 & ~isnan(sum(splitsInfo(:,2:end))) );
    splitsInfo = [splitsInfo(:,1) splitsInfo(:,1+fullIndx)];
    fullIndx = sum(splitsInfoSpace)~=0 & ~isnan(sum(splitsInfoSpace));
    splitsInfoSpace = splitsInfoSpace(:,fullIndx);
    splitsInfoSpaceExt = splitsInfoSpaceExt(:,fullIndx,:);
    fullIndx = sum(splitsInfoAmp)~=0 & ~isnan(sum(splitsInfoAmp));
    splitsInfoAmp = splitsInfoAmp(:,fullIndx);
end

%remove rows without merges or splits
if ~isempty(mergesInfo)
    filledRows = find(any(mergesInfo(:,2)~=0,2));
    mergesInfo = [filledRows mergesInfo(filledRows,:)];
    mergesInfoSpace = mergesInfoSpace(filledRows,:);
    mergesInfoSpaceExt = mergesInfoSpaceExt(filledRows,:,:);
    mergesInfoAmp = mergesInfoAmp(filledRows,:);
end
if ~isempty(splitsInfo)
    filledRows = find(any(splitsInfo(:,2)~=0,2));
    splitsInfo = [filledRows splitsInfo(filledRows,:)];
    splitsInfoSpace = splitsInfoSpace(filledRows,:);
    splitsInfoSpaceExt = splitsInfoSpaceExt(filledRows,:,:);
    splitsInfoAmp = splitsInfoAmp(filledRows,:);
end

%% Plotting

if plotRes
    switch probDim
        case 2
            plotMergeSplitPositions2D(tracks,mergesInfo,splitsInfo,...
                mergesInfoSpace,splitsInfoSpace,[],[],mergesInfoSpaceExt,splitsInfoSpaceExt)
    end
end

%% %%%%% ~~ the end ~~ %%%%%

