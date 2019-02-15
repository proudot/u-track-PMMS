function [msTimeMotionCharInfo] = calcMergeSplitTimesPartChar(tracks,...
    minTrackLen,probDim,diffAnalysisRes,removePotArtifacts)
%CALCMERGESPLITTIMESPARTCHAR calculates times between merges and splits and motion characteristics before and after merges and splits
%
%SYNOPSIS [msTimeMotionCharInfo] = calcMergeSplitTimesPartChar(tracks,...
%    minTrackLen,probDim,diffAnalysisRes,removePotArtifacts)
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
%OUTPUT msTimeInfo : Structure with field 'conf','brown' and 'linear' for
%                    confined, Brownian and linear tracks. Each field is a
%                    structure containing the fields:
%           .numTracks      : # of tracks in category.
%           .timeMerge2Split: 1-column vector storing the time interval
%                             between a merge and a consecutive split.
%           .timeSplit2MergeSelf: 1-column vector storing the time interval
%                             between a split and a consecutive merge 
%                             between the same two segments.
%           .timeSplit2MergeOther: 1-column vector storing the time interval
%                             between a split and a consecutive merge with
%                             different segments.
%           .timeMerge2End  : 1-column vector storing the time interval
%                             between a merge not followed by a split and
%                             the track's end.
%           .timeStart2Split: 1-column vector storing the time interval
%                             between a split not preceded by a merge and
%                             the track's start.
%           .charBeforeAfterMerge_timeMerge2Split: 4-column array with each
%                             row corresponding to the merge-to-split event
%                             stored in timeMerge2Split. The 4 columns show
%                             intensity before merging, intensity after
%                             merging, displacement before merging,
%                             displacement after merging. These values are
%                             calculated from 3 timepoints before and after
%                             each merge.
%           .charBeforeAfterSplit_timeStart2Split: 4-column array with each
%                             row corresponding to the start-to-split event
%                             stored in timeStart2Split. The 4 columns show
%                             intensity before splitting, intensity after
%                             splitting, displacement before splitting,
%                             displacement after splitting. These values are
%                             calculated from 3 timepoints before and after
%                             each split.
%           .charBeforeAfterMerge_timeMerge2End: 4-column array with each
%                             row corresponding to the merge-to-end event
%                             stored in timeMerge2End. The 4 columns show
%                             intensity before merging, intensity after
%                             merging, displacement before merging,
%                             displacement after merging. These values are
%                             calculated from 3 timepoints before and after
%                             each merge.
%                             
%
%Khuloud Jaqaman, November 2010

%% input

if nargin < 1 || isempty(tracks)
    disp('calcMergeSplitTimesPartChar: Missing input argument!');
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

%get number of tracks
numTracks = length(tracks);

%put tracks in matrix format for later use
[tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);
xCoordMat = tracksMat(:,1:8:end);
yCoordMat = tracksMat(:,2:8:end);
zCoordMat = tracksMat(:,3:8:end);
ampMat = tracksMat(:,4:8:end);

%get number of frames
numFrames = size(ampMat,2);

%% track types

%assign track types based on track segment types
%track type is taken as the type of the most dynamic segment
%linear > Brownian > confined
trackType = NaN(numTracks,1);
for iTrack = 1 : numTracks
    if any(diffAnalysisRes(iTrack).classification(:,1)==1)
        trackType(iTrack) = 3;
    else
        trackType(iTrack) = max(diffAnalysisRes(iTrack).classification(:,2));
    end
end

%store indices of tracks in the various categories
indxConf = find(trackType==1);
numTracksConf = length(indxConf);
indxBrown = find(trackType==2);
numTracksBrown = length(indxBrown);
indxLin = find(trackType==3);
numTracksLin = length(indxLin);

%% time from merges to splits and vice versa

for iType = 1 : 3

    switch iType
        case 1
            trackType = 'Conf';
        case 2
            trackType = 'Brown';
        case 3
            trackType = 'Lin';
    end

    %initialize some variables
    timeSplit2MergeSelf = [];
    timeSplit2MergeOther = [];
    timeMerge2Split = [];
    timeMerge2End = [];
    timeStart2Split = [];
    infoMerge2Split = [];
    infoMerge2End = [];
    infoStart2Split = [];
    eval(['indxTracks = indx' trackType ';']);

    %go over all tracks of this type ...
    for iTrack = indxTracks'

        %get track's sequence of events
        seqOfEvents = tracks(iTrack).seqOfEvents;

        %if requested, remove splits and merges that are most likely artifacts
        if removePotArtifacts
            seqOfEvents = removeSplitMergeArtifacts(seqOfEvents,1);
            seqOfEvents2 = removeSplitMergeArtifacts(seqOfEvents,0);
        end
        
        %% Splitting times
        
        %find indices where there are splits
        splitIndxGlob = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);

        %go over all splits and calculate split-to-merge times if the two
        %tracks merge back with each other
        splitIndxTmp = splitIndxGlob;
        takenMerge = [];
        pairsSplit2Merge = [];
        for iMS = 1 : length(splitIndxGlob)

            %get index of split
            iSplit = splitIndxGlob(iMS);

            %get time of split
            splitTime = seqOfEvents(iSplit,1);

            %get the indices of the splitting segments
            segmentIndx = seqOfEvents(iSplit,3:4);

            %check if these segments merge with each other at some point later
            iMerge = find(any(seqOfEvents(:,3:4)==segmentIndx(1),2) & ...
                any(seqOfEvents(:,3:4)==segmentIndx(2),2) & seqOfEvents(:,2)==2);

            if ~isempty(iMerge) %if they do merge with each other

                %calculate time of merge
                mergeTime = seqOfEvents(iMerge,1);

                %calculate split-to-merge time
                timeSplit2MergeTmp = mergeTime - splitTime;

                %store in global array
                timeSplit2MergeSelf = [timeSplit2MergeSelf; timeSplit2MergeTmp];

                %indicate index of participating split with a NaN
                splitIndxTmp(iMS) = NaN;

                %store indices of participating splits and merges
                takenMerge = [takenMerge; iMerge];
                pairsSplit2Merge = [pairsSplit2Merge; [iSplit iMerge]];

            end

        end %(for iMS = 1 : length(splitIndxGlob))

        %remove splits marked with a NaN
        splitIndxTmp = splitIndxTmp(~isnan(splitIndxTmp));

        %go over remaining splits and see if their segments merge
        %separately with other segments
        for iMS = 1 : length(splitIndxTmp)

            %get index of split
            iSplit = splitIndxTmp(iMS);

            %get time of split
            splitTime = seqOfEvents(iSplit,1);

            %get the indices of splitting segments
            segmentIndx = seqOfEvents(iSplit,3:4);

            %check whether segment 1 merges with a different segment
            iMerge = find(any(seqOfEvents(:,3:4)==segmentIndx(1),2) & ...
                ~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);
            iMerge = setdiff(iMerge,takenMerge);

            %get time of merge
            mergeTime = seqOfEvents(iMerge,1);

            %keep the earliest time of merge that comes after the time of
            %split
            iMerge = iMerge(mergeTime>splitTime);
            iMerge = min(iMerge);
            mergeTime = seqOfEvents(iMerge,1);

            %if segment 1 does merge with something else ...
            if ~isempty(iMerge)

                %calculate split-to-merge time
                timeSplit2MergeTmp = mergeTime - splitTime;

                %store in global array
                timeSplit2MergeOther = [timeSplit2MergeOther; timeSplit2MergeTmp];

                %store indices of participating splits and merges
                takenMerge = [takenMerge; iMerge];
                pairsSplit2Merge = [pairsSplit2Merge; [iSplit iMerge]];

            end

            %check whether segment 2 merges with a different segment
            iMerge = find(any(seqOfEvents(:,3:4)==segmentIndx(2),2) & ...
                ~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);
            iMerge = setdiff(iMerge,takenMerge);

            %get time of merge
            mergeTime = seqOfEvents(iMerge,1);

            %keep the earliest time of merge that comes after the time of
            %split
            iMerge = iMerge(mergeTime>splitTime);
            iMerge = min(iMerge);
            mergeTime = seqOfEvents(iMerge,1);

            %if segment 2 does merge with something else ...
            if ~isempty(iMerge)

                %calculate split-to-merge time
                timeSplit2MergeTmp = mergeTime - splitTime;

                %store in global array
                timeSplit2MergeOther = [timeSplit2MergeOther; timeSplit2MergeTmp];

                %store indices of participating splits and merges
                takenMerge = [takenMerge; iMerge];
                pairsSplit2Merge = [pairsSplit2Merge; [iSplit iMerge]];

            end

        end %(for iMS = 1 : length(splitIndxTmp))
        
        %% Merging times

        %find indices where there are merges
        mergeIndxGlob = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);

        %go over all merges and see if there are consequent splits
        takenSplit = [];
        pairsMerge2Split = [];
        for iMS = 1 : length(mergeIndxGlob)

            %get index of merge
            iMerge = mergeIndxGlob(iMS);

            %get time of merge
            mergeTime = seqOfEvents(iMerge,1);

            %get indices of merging segments
            segmentIndx = seqOfEvents(iMerge,3:4);

            %check wether the continuing segment splits again
            iSplit = find(any(seqOfEvents(:,3:4)==segmentIndx(2),2) & ...
                ~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
            iSplit = setdiff(iSplit,takenSplit);

            %get times of split
            splitTime = seqOfEvents(iSplit,1);

            %keep the earliest time of split that comes after the time of
            %merge
            iSplit = iSplit(splitTime>mergeTime);
            iSplit = min(iSplit);
            splitTime = seqOfEvents(iSplit,1);

            %if there is a consequent split ...
            if ~isempty(iSplit)

                %calculate the merge-to-split time
                timeMerge2SplitTmp = splitTime - mergeTime;

                %store in global array
                timeMerge2Split = [timeMerge2Split; timeMerge2SplitTmp];

                %store indices of participating merges and splits
                takenSplit = [takenSplit; iSplit];
                pairsMerge2Split = [pairsMerge2Split; [iMerge iSplit]];
                
                %to calculate particle characteristics before and after the
                %merge ...
                
                %get indices of participating segments in the global segment matrix
                segmentsMS = seqOfEvents2(iMerge,3:4) + trackStartRow(iTrack) - 1;
                
                %calculate some indices taking into account movie start and end
                %times
                msTimeMinus3 = max(mergeTime-3,1);
                msTimeMinus1 = max(mergeTime-1,1);
                msTimePlus2  = min(mergeTime+2,numFrames);
                
                %calculate the intensity characteristics before and after
                %the merge
                intBefore = ampMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                intBefore = nanmean(intBefore(:));
                intAfter  = ampMat(segmentsMS,mergeTime:msTimePlus2);
                intAfter  = nanmean(intAfter(:));
                
                %calculate the displacement characteristics before the
                %merge
                xBefore = xCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                yBefore = yCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                zBefore = zCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                dispBefore = sqrt(diff(xBefore,[],2).^2 + diff(yBefore,[],2).^2 + ...
                    diff(zBefore,[],2).^2);
                dispBefore = nanmean(dispBefore(:));
                
                %calculate the displacement characteristics after the merge
                xAfter  = xCoordMat(segmentsMS,mergeTime:msTimePlus2);
                yAfter  = yCoordMat(segmentsMS,mergeTime:msTimePlus2);
                zAfter  = zCoordMat(segmentsMS,mergeTime:msTimePlus2);
                dispAfter = sqrt(diff(xAfter,[],2).^2 + diff(yAfter,[],2).^2 + ...
                    diff(zAfter,[],2).^2);
                dispAfter = nanmean(dispAfter(:));
                
                %save this information
                infoMerge2Split = [infoMerge2Split; [intBefore intAfter dispBefore dispAfter]];
                
            end
            
        end %(for iMS = 1 : length(mergeIndxGlob))
        
        %go over all splits and calculate start-to-split times if appropriate
        for iMS = 1 : length(splitIndxGlob)
            
            %get index of split
            iSplit = splitIndxGlob(iMS);

            %check if split is preceded by a merge
            if ~isempty(pairsMerge2Split)
                splitAfterMerge = length(find(pairsMerge2Split(:,2)==iSplit));
            else
                splitAfterMerge = 0;
            end

            %if split is not preceded by a merge, calculate start-to-split time
            if splitAfterMerge == 0

                %get time of split
                splitTime = seqOfEvents(iSplit,1);

                %get the index of the segment that got split out of
                segmentIndx = seqOfEvents(iSplit,4);

                %get the start time of this segment
                startTime = seqOfEvents(seqOfEvents(:,2)==1&...
                    seqOfEvents(:,3)==segmentIndx&isnan(seqOfEvents(:,4)),1);

                %calculate the start-to-split time
                timeStart2SplitTmp = splitTime - startTime;
                
                if ~isempty(timeStart2SplitTmp)
                    
                    %store in global array
                    timeStart2Split = [timeStart2Split; timeStart2SplitTmp];
                    
                    %to calculate particle characteristics before and after the
                    %split ...
                    
                    %get indices of participating segments in the global segment matrix
                    segmentsMS = seqOfEvents2(iSplit,3:4) + trackStartRow(iTrack) - 1;
                    
                    %calculate some indices taking into account movie start and end
                    %times
                    msTimeMinus3 = max(splitTime-3,1);
                    msTimeMinus1 = max(splitTime-1,1);
                    msTimePlus2  = min(splitTime+2,numFrames);
                    
                    %calculate the intensity characteristics before and after
                    %the split
                    intBefore = ampMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    intBefore = nanmean(intBefore(:));
                    intAfter  = ampMat(segmentsMS,splitTime:msTimePlus2);
                    intAfter  = nanmean(intAfter(:));
                    
                    %calculate the displacement characteristics before the
                    %split
                    xBefore = xCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    yBefore = yCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    zBefore = zCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    dispBefore = sqrt(diff(xBefore,[],2).^2 + diff(yBefore,[],2).^2 + ...
                        diff(zBefore,[],2).^2);
                    dispBefore = nanmean(dispBefore(:));
                    
                    %calculate the displacement characteristics after the split
                    xAfter  = xCoordMat(segmentsMS,splitTime:msTimePlus2);
                    yAfter  = yCoordMat(segmentsMS,splitTime:msTimePlus2);
                    zAfter  = zCoordMat(segmentsMS,splitTime:msTimePlus2);
                    dispAfter = sqrt(diff(xAfter,[],2).^2 + diff(yAfter,[],2).^2 + ...
                        diff(zAfter,[],2).^2);
                    dispAfter = nanmean(dispAfter(:));
                    
                    %save this information
                    infoStart2Split = [infoStart2Split; [intBefore intAfter dispBefore dispAfter]];
                    
                end

            end

        end %(for iMS = 1 : length(splitIndxGlob))

        %go over all merges and calculate merge-to-end times if appropriate
        for iMS = 1 : length(mergeIndxGlob)

            %get index of merge
            iMerge = mergeIndxGlob(iMS);

            %check if merge is followed by a split
            if ~isempty(pairsMerge2Split)
                splitAfterMerge = length(find(pairsMerge2Split(:,1)==iMerge));
            else
                splitAfterMerge = 0;
            end

            %if merge is not followed by a split, calculate merge-to-end time
            if splitAfterMerge == 0

                %get time of merge
                mergeTime = seqOfEvents(iMerge,1);

                %get the index of the segment that got merged with
                segmentIndx = seqOfEvents(iMerge,4);

                %get the end time of this segment
                endTime = seqOfEvents(seqOfEvents(:,2)==2&...
                    seqOfEvents(:,3)==segmentIndx&isnan(seqOfEvents(:,4)),1);

                %calculate the merge-to-end time
                timeMerge2EndTmp = endTime - mergeTime;
                
                if ~isempty(timeMerge2EndTmp)
                    
                    %store in global array
                    timeMerge2End = [timeMerge2End; timeMerge2EndTmp];
                    
                    %to calculate particle characteristics before and after the
                    %merge ...
                    
                    %get indices of participating segments in the global segment matrix
                    segmentsMS = seqOfEvents2(iMerge,3:4) + trackStartRow(iTrack) - 1;
                    
                    %calculate some indices taking into account movie start and end
                    %times
                    msTimeMinus3 = max(mergeTime-3,1);
                    msTimeMinus1 = max(mergeTime-1,1);
                    msTimePlus2 = min(mergeTime+2,numFrames);
                    
                    %calculate the intensity characteristics before and after
                    %the merge
                    intBefore = ampMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    intBefore = nanmean(intBefore(:));
                    intAfter  = ampMat(segmentsMS,mergeTime:msTimePlus2);
                    intAfter  = nanmean(intAfter(:));
                    
                    %calculate the displacement characteristics before the
                    %merge
                    xBefore = xCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    yBefore = yCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    zBefore = zCoordMat(segmentsMS,msTimeMinus3:msTimeMinus1);
                    dispBefore = sqrt(diff(xBefore,[],2).^2 + diff(yBefore,[],2).^2 + ...
                        diff(zBefore,[],2).^2);
                    dispBefore = nanmean(dispBefore(:));
                    
                    %calculate the displacement characteristics after the merge
                    xAfter  = xCoordMat(segmentsMS,mergeTime:msTimePlus2);
                    yAfter  = yCoordMat(segmentsMS,mergeTime:msTimePlus2);
                    zAfter  = zCoordMat(segmentsMS,mergeTime:msTimePlus2);
                    dispAfter = sqrt(diff(xAfter,[],2).^2 + diff(yAfter,[],2).^2 + ...
                        diff(zAfter,[],2).^2);
                    dispAfter = nanmean(dispAfter(:));
                    
                    %save this information
                    infoMerge2End = [infoMerge2End; [intBefore intAfter dispBefore dispAfter]];
                    
                end
                
            end

        end %(for iMS = 1 : length(mergeIndxGlob))

    end %(for iTrack = indxTracks')

    %store track numbers and distributions
    %     eval(['numTracks' trackType ' = [numTrackMS numTrackOnlyM numTrackOnlyS numTrackNoMS];'])
    eval(['timeMerge2Split' trackType ' = timeMerge2Split;'])
    eval(['timeSplit2MergeSelf' trackType ' = timeSplit2MergeSelf;'])
    eval(['timeSplit2MergeOther' trackType ' = timeSplit2MergeOther;'])
    eval(['timeMerge2End' trackType ' = timeMerge2End;'])
    eval(['timeStart2Split' trackType ' = timeStart2Split;'])
    eval(['charBeforeAfterMerge_timeMerge2Split' trackType ' = infoMerge2Split;'])
    eval(['charBeforeAfterSplit_timeStart2Split' trackType ' = infoStart2Split;'])
    eval(['charBeforeAfterMerge_timeMerge2End' trackType ' = infoMerge2End;'])

end

%% output

msTimeInfo.linear.numTracks = numTracksLin;
msTimeInfo.linear.timeMerge2Split = timeMerge2SplitLin;
msTimeInfo.linear.timeSplit2MergeSelf = timeSplit2MergeSelfLin;
msTimeInfo.linear.timeSplit2MergeOther = timeSplit2MergeOtherLin;
msTimeInfo.linear.timeStart2Split = timeStart2SplitLin;
msTimeInfo.linear.timeMerge2End = timeMerge2EndLin;
msTimeInfo.linear.charBeforeAfterMerge_timeMerge2Split = ...
    charBeforeAfterMerge_timeMerge2SplitLin;
msTimeInfo.linear.charBeforeAfterSplit_timeStart2Split = ...
    charBeforeAfterSplit_timeStart2SplitLin;
msTimeInfo.linear.charBeforeAfterMerge_timeMerge2End = ...
    charBeforeAfterMerge_timeMerge2EndLin;

msTimeInfo.brown.numTracks = numTracksBrown;
msTimeInfo.brown.timeMerge2Split = timeMerge2SplitBrown;
msTimeInfo.brown.timeSplit2MergeSelf = timeSplit2MergeSelfBrown;
msTimeInfo.brown.timeSplit2MergeOther = timeSplit2MergeOtherBrown;
msTimeInfo.brown.timeStart2Split = timeStart2SplitBrown;
msTimeInfo.brown.timeMerge2End = timeMerge2EndBrown;
msTimeInfo.brown.charBeforeAfterMerge_timeMerge2Split = ...
    charBeforeAfterMerge_timeMerge2SplitBrown;
msTimeInfo.brown.charBeforeAfterSplit_timeStart2Split = ...
    charBeforeAfterSplit_timeStart2SplitBrown;
msTimeInfo.brown.charBeforeAfterMerge_timeMerge2End = ...
    charBeforeAfterMerge_timeMerge2EndBrown;

msTimeInfo.conf.numTracks = numTracksConf;
msTimeInfo.conf.timeMerge2Split = timeMerge2SplitConf;
msTimeInfo.conf.timeSplit2MergeSelf = timeSplit2MergeSelfConf;
msTimeInfo.conf.timeSplit2MergeOther = timeSplit2MergeOtherConf;
msTimeInfo.conf.timeStart2Split = timeStart2SplitConf;
msTimeInfo.conf.timeMerge2End = timeMerge2EndConf;
msTimeInfo.conf.charBeforeAfterMerge_timeMerge2Split = ...
    charBeforeAfterMerge_timeMerge2SplitConf;
msTimeInfo.conf.charBeforeAfterSplit_timeStart2Split = ...
    charBeforeAfterSplit_timeStart2SplitConf;
msTimeInfo.conf.charBeforeAfterMerge_timeMerge2End = ...
    charBeforeAfterMerge_timeMerge2EndConf;

msTimeMotionCharInfo = msTimeInfo;

%% ~~~ the end ~~~


% % %     %give splits a negative time and thus keep only time vector
% % %     msTime(msTime(:,2)==1,1) = -msTime(msTime(:,2)==1,1);
% % %     msTime = msTime(:,1);
% % %
% % %         %find the start time and end time of the compound track
% % %         startTime = min(seqOfEvents(:,1));
% % %         endTime = max(seqOfEvents(:,1));
% % %
% % %         %take action based on whether there are merges and splits
% % %         if isempty(msTime) %if there are no merges and no splits
% % %
% % %             %add one to counter of tracks without merges and splits
% % %             numTrackNoMS = numTrackNoMS + 1;
% % %
% % %         elseif all(msTime > 0) %if there are merges but no splits
% % %
% % %             %add one to counter of tracks with only merges
% % %             numTrackOnlyM = numTrackOnlyM + 1;
% % %
% % %             %consider time intervel from time of merging to end as "merging
% % %             %time"
% % %             timeMerge2End = [timeMerge2End; endTime-msTime];
% % %
% % %         elseif all(msTime < 0) %if there are splits but no merges
% % %
% % %             %add one to counter of tracks with only splits
% % %             numTrackOnlyS = numTrackOnlyS + 1;
% % %
% % %             %consider time intervel from start to time of splitting as
% % %             %"merging time"
% % %             timeStart2Split = [timeStart2Split; abs(msTime)-startTime];
% % %
% % %         else %if there are both merges and splits
% % %
% % %             %add one to counter of tracks with both merges and
% % %             %splits
% % %             numTrackMS = numTrackMS + 1;
% % %
% % %             %get number of merges and splits to consider
% % %             numMS = length(msTime);
% % %
% % %             indxMS = 1;
% % %             while indxMS < numMS
% % %
% % %                 %get index of initial event and its type
% % %                 typeEvent1 = sign(msTime(indxMS));
% % %                 indxEvent1 = indxMS;
% % %
% % %                 %keep increasing indxMS until you find an event different
% % %                 %from initial event
% % %                 typeEvent2 = typeEvent1;
% % %                 while typeEvent2 == typeEvent1 && indxMS < numMS
% % %                     indxMS = indxMS + 1;
% % %                     typeEvent2 = sign(msTime(indxMS));
% % %                 end
% % %                 indxEvent2 = indxMS;
% % %
% % %                 %keep increasing indxMS until you find an event again
% % %                 %similar to initial event
% % %                 typeEvent3 = typeEvent2;
% % %                 while typeEvent3 ~= typeEvent1 && indxMS < numMS
% % %                     indxMS = indxMS + 1;
% % %                     typeEvent3 = sign(msTime(indxMS));
% % %                 end
% % %                 if typeEvent3 ~= typeEvent2
% % %                     indxMS = indxMS - 1;
% % %                 end
% % %                 indxEvent3 = indxMS;
% % %
% % %                 if typeEvent2 ~= typeEvent1
% % %
% % %                     %calculate total number of combinations for calculating
% % %                     %time between merges and splits (or vice versa)
% % %                     numCombination = (indxEvent2-indxEvent1)*(indxEvent3-indxEvent2+1);
% % %
% % %                     %go over all combinations and calculate time
% % %                     indxComb = 0;
% % %                     timeBetweenMS = zeros(numCombination,1);
% % %                     for indx1 = 1 : indxEvent2-indxEvent1
% % %                         for indx2 = 1 : indxEvent3-indxEvent2+1
% % %                             indxComb = indxComb + 1;
% % %                             timeBetweenMS(indxComb) = abs(msTime(indx2+indxEvent2-1)) - ...
% % %                                 abs(msTime(indx1+indxEvent1-1));
% % %                         end
% % %                     end
% % %
% % %                     %store the times and their weights based on whether we
% % %                     %looked at a merge to split or a split to merge
% % %                     if typeEvent1 > 0
% % %                         timeMerge2Split = [timeMerge2Split; [timeBetweenMS ...
% % %                             (1/numCombination)*ones(numCombination,1)]];
% % %                     else
% % %                         timeSplit2Merge = [timeSplit2Merge; [timeBetweenMS ...
% % %                             (1/numCombination)*ones(numCombination,1)]];
% % %                     end
% % %
% % %                 end
% % %
% % %                 %update indxMS to look at next merge-to-split or
% % %                 %split-to-merge event
% % %                 indxMS = indxEvent2;
% % %
% % %             end %(while indxMS <= numMS - 1)
% % %
% % %             %if initial event is a split, consider time from beginning to
% % %             %time of split as merging time
% % %             indxMS = 1;
% % %             typeEvent1 = sign(msTime(indxMS));
% % %             if typeEvent1 < 0
% % %
% % %                 %keep increasing indxMS until you find an event different
% % %                 %from initial event
% % %                 typeEvent2 = typeEvent1;
% % %                 while typeEvent2 == typeEvent1
% % %                     indxMS = indxMS + 1;
% % %                     typeEvent2 = sign(msTime(indxMS));
% % %                 end
% % %
% % %                 %add time intervals to "merging time"
% % %                 timeStart2Split = [timeStart2Split; abs(msTime(1:indxMS-1))-startTime];
% % %
% % %             end
% % %
% % %             %if final event is a merge, consider time from time of merge to
% % %             %end as merging time
% % %             indxMS = numMS;
% % %             typeEvent1 = sign(msTime(indxMS));
% % %             if typeEvent1 > 0
% % %
% % %                 %keep decreasing indxMS until you find an event different
% % %                 %from final event
% % %                 typeEvent2 = typeEvent1;
% % %                 while typeEvent2 == typeEvent1
% % %                     indxMS = indxMS - 1;
% % %                     typeEvent2 = sign(msTime(indxMS));
% % %                 end
% % %
% % %                 %add time intervals to "merging time"
% % %                 timeMerge2End = [timeMerge2End; endTime-msTime(indxMS+1:numMS)];
% % %
% % %             end
% % %
% % %         end %(if isempty(msTime) ... elseif ...)
% % %
% % %     end %(for iTrack = 1 : indxTracks')
% % %
