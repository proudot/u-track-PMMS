function [msTimeInfo] = calcMergeSplitTimes(tracks,minTrackLen,...
    diffAnalysisRes,removePotArtifacts)
%CALCMERGESPLITTIMES calculates times between merges and splits
%
%SYNOPSIS [msTimeInfo] = calcMergeSplitTimes(tracks,minTrackLen,...
%    diffAnalysisRes,removePotArtifacts)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffAnalysis1).
%                    Optional. If not input, all tracks will be analyzed
%                    together.
%       removePotArtifacts: 1 to remove potentially artifactual merges and
%                    splits, resulting for instance from detection
%                    artifact, 0 otherwise. 
%                    Optional. Default: 1.
%
%OUTPUT msTimeInfo : Structure with fields 'conf','brown' and 'linear' for
%                    confined, Brownian and linear tracks OR field 'all' if
%                    diffAnalysisRes was not input. Each field is a
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
%
%Khuloud Jaqaman, April 2008

%% input

if nargin < 1 || isempty(tracks)
    disp('calcStatsMS: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(diffAnalysisRes)
    diffAnalysisRes = [];
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
if ~isempty(diffAnalysisRes)
    diffAnalysisRes = diffAnalysisRes(indx);
end

%get number of tracks and number of frames
numTracks = length(tracks);

%% track types

if ~isempty(diffAnalysisRes)
    
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
    
    %determine number of track types
    numType = 3;
    
else
    
    indxAll = (1:numTracks)';
    numTracksAll = length(indxAll);
    numType = 1;
    
end

%% time from merges to splits and vice versa

for iType = 1 : numType

    if numType == 3
        switch iType
            case 1
                trackType = 'Conf';
            case 2
                trackType = 'Brown';
            case 3
                trackType = 'Lin';
        end
    else
        trackType = 'All';
    end

    %initialize some variables
    timeSplit2MergeSelf = [];
    timeSplit2MergeOther = [];
    timeMerge2Split = [];
    timeMerge2End = [];
    timeStart2Split = [];
    eval(['indxTracks = indx' trackType ';']);

    %go over all tracks of this type ...
    for iTrack = indxTracks'

        %get track's sequence of events
        seqOfEvents = tracks(iTrack).seqOfEvents;

        %if requested, remove splits and merges that are most likely artifacts
        if removePotArtifacts
            seqOfEvents = removeSplitMergeArtifacts(seqOfEvents,1);
        end
        
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

                %store in global array
                timeStart2Split = [timeStart2Split; timeStart2SplitTmp];

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

                %store in global array
                timeMerge2End = [timeMerge2End; timeMerge2EndTmp];

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

end

%% output

if ~isempty(diffAnalysisRes)
    msTimeInfo.linear.numTracks = numTracksLin;
    msTimeInfo.linear.timeMerge2Split = timeMerge2SplitLin;
    msTimeInfo.linear.timeSplit2MergeSelf = timeSplit2MergeSelfLin;
    msTimeInfo.linear.timeSplit2MergeOther = timeSplit2MergeOtherLin;
    msTimeInfo.linear.timeStart2Split = timeStart2SplitLin;
    msTimeInfo.linear.timeMerge2End = timeMerge2EndLin;
    msTimeInfo.brown.numTracks = numTracksBrown;
    msTimeInfo.brown.timeMerge2Split = timeMerge2SplitBrown;
    msTimeInfo.brown.timeSplit2MergeSelf = timeSplit2MergeSelfBrown;
    msTimeInfo.brown.timeSplit2MergeOther = timeSplit2MergeOtherBrown;
    msTimeInfo.brown.timeStart2Split = timeStart2SplitBrown;
    msTimeInfo.brown.timeMerge2End = timeMerge2EndBrown;
    msTimeInfo.conf.numTracks = numTracksConf;
    msTimeInfo.conf.timeMerge2Split = timeMerge2SplitConf;
    msTimeInfo.conf.timeSplit2MergeSelf = timeSplit2MergeSelfConf;
    msTimeInfo.conf.timeSplit2MergeOther = timeSplit2MergeOtherConf;
    msTimeInfo.conf.timeStart2Split = timeStart2SplitConf;
    msTimeInfo.conf.timeMerge2End = timeMerge2EndConf;
else
    msTimeInfo.all.numTracks = numTracksAll;
    msTimeInfo.all.timeMerge2Split = timeMerge2SplitAll;
    msTimeInfo.all.timeSplit2MergeSelf = timeSplit2MergeSelfAll;
    msTimeInfo.all.timeSplit2MergeOther = timeSplit2MergeOtherAll;
    msTimeInfo.all.timeStart2Split = timeStart2SplitAll;
    msTimeInfo.all.timeMerge2End = timeMerge2EndAll;
end

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
