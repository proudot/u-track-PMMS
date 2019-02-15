function [statsPerTrack] = scoreTrackLinksGapsMSLft(tracksFinal,tracksSim)
%SCORETRACKLINKSGAPSMSLFT evaluates links, gaps, merges/splits and lifetime per track
%
%SYNOPSIS [statsPerTrack] = scoreTrackLinksGapsMSLft(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : Simulated tracks as obtained from
%                     simulateMimickCD36_MS (structure format).
%
%OUTPUT statsPerTrack: Number of tracks - by - 5 array. Row i
%                      corresponds to track i. The columns show:
%                      (1) Number of wrong links.
%                      (2) Number of wrongly closed gaps.
%                      (3) Number of wrong merges/splits.
%                      (4) Lifetime of track.
%                      (5) Lifetime of corresponding ground truth track,
%                          when there is a ground truth track that starts
%                          at the same place and in the same frame. If
%                          there is no corresponding ground truth track,
%                          -15 will be shown instead.
%                      In (1)-(3), false positives implicitly contribute to
%                      wrong links, gaps and merges/split.
%
%REMARKS I'm not sure how well the lifetime calculation works when there
%        are merges and splits.
%
%Khuloud Jaqaman, October 2007

%% process input variables

%extract track information out of tracksSim
[xyCoordAll0,xyCoordM0,xyCoordS0] = trackInfoFromStruct(tracksSim);

%get ground truth track start, end and life times
trackSEL0 = getTrackSEL(tracksSim);

%get coordinates of each track at its start
trackCoordStart0 = NaN(length(tracksSim),2);
for iTrack = 1 : length(tracksSim)
    trackCoordStart0(iTrack,:) = tracksSim(iTrack).tracksCoordAmpCG(...
        tracksSim(iTrack).seqOfEvents(1,3),1:2);
end

if isstruct(tracksFinal)

    %extract track information out of tracksFinal
    [xyCoordAll1,xyCoordM1,xyCoordS1,trackStartRow1] = trackInfoFromStruct(tracksFinal);
    trackStartRow1(end+1) = size(xyCoordAll1,1) + 1;

    %get number of tracks
    numTracks = length(tracksFinal);

    %get track start, end and life times
    trackSEL1 = getTrackSEL(tracksFinal);

    %get coordinates of each track at its start
    trackCoordStart1 = NaN(numTracks,2);
    for iTrack = 1 : numTracks
        trackCoordStart1(iTrack,:) = tracksFinal(iTrack).tracksCoordAmpCG(...
            tracksFinal(iTrack).seqOfEvents(1,3),1:2);
    end

else

    %get number of tracks
    numTracks = size(tracksFinal,1);

    %directly extract track information out of tracksFinal
    xyCoordAll1 = zeros(numTracks,size(tracksFinal,2)/4);
    xyCoordAll1(:,1:2:end) = tracksFinal(:,1:8:end);
    xyCoordAll1(:,2:2:end) = tracksFinal(:,2:8:end);
    xyCoordAll1(isnan(xyCoordAll1)) = 0;
    trackStartRow1 = (1:numTrack+1)';
    xyCoordS1_all
    %indicate that there are no merges and splits
    xyCoordM1 = [];
    xyCoordS1 = [];

    %get track start, end and life times
    trackSEL1 = getTrackSEL(tracksFinal);

    %get coordinates of each track at its start
    trackCoordStart1 = NaN(numTracks,2);
    for iTrack = 1 : numTracks
        trackCoordStart1(iTrack,:) = tracksFinal(iTrack,...
            (trackSEL1(iTrack,1)-1)*8+1:(trackSEL1(iTrack,1)-1)*8+2);
    end

end

%% check links, gaps, merges/splits and lifetimes

%initialize output variable
statsPerTrack = NaN(numTracks,5);

%go over all tracks ...
for iTrack = 1 : numTracks


    %% first, evaluate this track's links and gaps

    %initialize to zero number of wrong links and number of wrong gaps
    numLinksWrong = 0;
    numGapsWrong = 0;

    %go over all segments (rows) making this track
    for iRow = trackStartRow1(iTrack) : trackStartRow1(iTrack+1)-1

        %find frames where this track exists
        framesExist = find(xyCoordAll1(iRow,1:2:end));

        %go over these frames ...
        for iIndx = 1 : length(framesExist) - 1

            %get current frame and next frame
            iFrame1 = framesExist(iIndx);
            iFrame2 = framesExist(iIndx+1);

            %get coordinates in these two frames
            xyCoord12 = [xyCoordAll1(iRow,2*iFrame1-1:2*iFrame1); ...
                xyCoordAll1(iRow,2*iFrame2-1:2*iFrame2)];

            %find location of coordinates in first frame in ground truth
            rowsGT = find( xyCoordAll0(:,iFrame1*2-1)==xyCoord12(1,1) & ...
                xyCoordAll0(:,iFrame1*2)==xyCoord12(1,2) );

            %if coordinates in first frame exist in ground truth ...
            if ~isempty(rowsGT)

                %get corresponding coordinates in second frame in ground truth
                xyCoord02 = [xyCoordAll0(rowsGT,iFrame2*2-1) xyCoordAll0(rowsGT,iFrame2*2)];

                %check whether they are the same as the coordinates in the second
                %frame in the tracking results
                commonCoord = intersect(xyCoord02,xyCoord12(2,:),'rows');

                %if the coordinates in the second frame are different between the
                %tracking results and the ground truth, then link/closed gap is wrong
                if isempty(commonCoord)
                    if (iFrame2 - iFrame1) == 1
                        numLinksWrong = numLinksWrong + 1;
                    else
                        numGapsWrong = numGapsWrong + 1;
                    end
                end

            else %if coordinates do not exist in g.t., i.e. they represent a false positive

                %this link/closed gap is also wrong
                if (iFrame2 - iFrame1) == 1
                    numLinksWrong = numLinksWrong + 1;
                else
                    numGapsWrong = numGapsWrong + 1;
                end

            end

        end %(for iIndx = 1 : length(framesExist) - 1)

    end %(for iRow = trackStartRow1(iTrack) : trackStartRow1(iTrack+1)-1)


    %% second, evaluate this track's merges and splits

    %find merging and splitting events in this track
    if isempty(xyCoordM1)
        indxMerge = [];
    else
        indxMerge = find(xyCoordM1(:,1) == iTrack);
    end
    if isempty(xyCoordS1)
        indxSplit = [];
    else
        indxSplit = find(xyCoordS1(:,1) == iTrack);
    end

    %check whether there are merges
    if isempty(indxMerge)

        %if there are no merges, this calculation is irrelevant
        numMergesWrong = NaN;

    else %if there are merges

        %go over merges
        numMergesWrong = 0;
        for iMerge = indxMerge'

            %get time of merge
            timeMerge = xyCoordM1(iMerge,2);

            %get coordinates before and after merge
            xyCoordMAft1  = xyCoordM1(iMerge,3:4);
            xyCoordMBef11 = xyCoordM1(iMerge,5:6);
            xyCoordMBef12 = xyCoordM1(iMerge,7:8);

            %check if this merge exists in the ground truth
            [dummy,dummy,indxGT] = intersect([timeMerge xyCoordMAft1],xyCoordM0(:,2:4),'rows');

            %if merge does not exist in ground truth
            if isempty(indxGT)

                %add 1 to number of wrong merges
                numMergesWrong = numMergesWrong + 1;

            else %if merge exists in ground truth

                %get coordinates of merging features in ground truth
                xyCoordMBef01 = xyCoordM0(indxGT,5:6);
                xyCoordMBef02 = xyCoordM0(indxGT,7:8);

                %find how many tracking and ground truth coordinates are equivalent
                numEquiv = size(intersect([xyCoordMBef11; xyCoordMBef12],...
                    [xyCoordMBef01; xyCoordMBef02],'rows'),1);

                %update number of wrong merges based on numEquiv
                switch numEquiv

                    case 2 %if both coordinates have equivalents

                        %this is a correct merge

                    case 1 %if only one coordinate has an equivalent

                        %add 1/2 to number of wrong merges
                        numMergesWrong = numMergesWrong + 0.5;

                    case 0 %if neither coordinate has an equivalent

                        %add 1 to number of wrong merges
                        numMergesWrong = numMergesWrong + 1;

                end %(switch numEquiv)

            end %(if isempty(indxGT) ... else ...)

        end %(for iMerge = indxMerge')

    end %(if isempty(indxMerge) ... else ...)

    %check whether there are splits
    if isempty(indxSplit)

        %if there are no splits, this calculation is irrelevant
        numSplitsWrong = NaN;

    else %if there are splits

        %go over splits
        numSplitsWrong = 0;
        for iSplit = indxSplit'

            %get time of split
            timeSplit = xyCoordS1(iSplit,2);

            %get coordinates before and after split
            xyCoordSBef1  = xyCoordS1(iSplit,3:4);
            xyCoordSAft11 = xyCoordS1(iSplit,5:6);
            xyCoordSAft12 = xyCoordS1(iSplit,7:8);

            %check if this split exists in the ground truth
            [dummy,dummy,indxGT] = intersect([timeSplit xyCoordSBef1],xyCoordS0(:,2:4),'rows');

            %if split does not exist in ground truth
            if isempty(indxGT)

                %add 1 to number of wrong splits
                numSplitsWrong = numSplitsWrong + 1;

            else %if split exists in ground truth

                %get coordinates of splitting features in ground truth
                xyCoordSAft01 = xyCoordS0(indxGT,5:6);
                xyCoordSAft02 = xyCoordS0(indxGT,7:8);

                %find how many tracking and ground truth coordinates are equivalent
                numEquiv = size(intersect([xyCoordSAft11; xyCoordSAft12],...
                    [xyCoordSAft01; xyCoordSAft02],'rows'),1);

                %update number of wrong splits based on numEquiv
                switch numEquiv

                    case 2 %if both coordinates have equivalents

                        %this is a correct split

                    case 1 %if only one coordinate has an equivalent

                        %add 1/2 to number of wrong splits
                        numSplitsWrong = numSplitsWrong + 0.5;

                    case 0 %if neither coordinate has an equivalent

                        %add 1 to number of wrong splits
                        numSplitsWrong = numSplitsWrong + 1;

                end %(switch numEquiv)

            end %(if isempty(indxGT) ... else ...)

        end %(for iSplit = indxSplit')

    end %(if isempty(indxSplit) ... else ...)

    if isnan(numMergesWrong) && isnan(numSplitsWrong)
        numMSWrong = NaN;
    else
        numMSWrong = nansum([numMergesWrong numSplitsWrong]);
    end


    %% third, if the start of this track corresponds to the start of
    %% a track in the ground truth, compare track's lifetime to that 
    %% of the ground truth track

    %get frame and coordinates at start of this track
    iFrame1 = trackSEL1(iTrack,1);
    xCoordStart1 = trackCoordStart1(iTrack,1);
    yCoordStart1 = trackCoordStart1(iTrack,2);

    %get track's lifetime
    trackLft = trackSEL1(iTrack,3);

    %look for a ground truth track that starts at the same frame with the
    %same coordinates
    gtTrack = find(trackSEL0(:,1)==iFrame1 & trackCoordStart0(:,1)==xCoordStart1 & ...
        trackCoordStart0(:,2)==yCoordStart1);

    %if there is no corresponding ground truth track ...
    if isempty(gtTrack)

        %assign -15 to ground truth lifetime
        groundTruthLft = -15;

    else %if there is a corresponding ground truth track

        %get its lifetime
        groundTruthLft = trackSEL0(gtTrack,3);

    end

    %% store this track's information

    statsPerTrack(iTrack,:) = [numLinksWrong numGapsWrong numMSWrong ...
        trackLft groundTruthLft];

end %(for iTrack = 1 : numTracks)


%% Subfunction 1

function [xyCoordAll,xyCoordM,xyCoordS,trackStartRow] = trackInfoFromStruct(tracks)

%get number of tracks
numTracks = length(tracks);

%get number of frames
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%initialize output variables
xyCoordAll = [];
xyCoordM = [];
xyCoordS = [];
trackStartRow = zeros(numTracks,1);

%go over all tracks ...
for i = 1 : numTracks

    %get sequence of events of track
    seqOfEvents = tracks(i).seqOfEvents;

    %get start and end times of track
    startTime = seqOfEvents(1,1);
    endTime   = seqOfEvents(end,1);

    %extract track coordinates from structure
    tracksCoordAmpCG = tracks(i).tracksCoordAmpCG;
    xyCoordTmp = zeros(size(tracksCoordAmpCG,1),numFrames*2);
    xyCoordTmp(:,2*(startTime-1)+1:2:2*(endTime-1)+1) = tracksCoordAmpCG(:,1:8:end);
    xyCoordTmp(:,2*(startTime-1)+2:2:2*(endTime-1)+2) = tracksCoordAmpCG(:,2:8:end);

    %determine row where this compound track starts in xyCoordAll
    trackStartRow(i) = size(xyCoordAll,1) + 1;

    %get merge events
    mergeIndx = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));

    %go over all merges
    for iMerge = mergeIndx'

        %get merge time
        timeMerge = seqOfEvents(iMerge,1);

        %store coordinates belonging to this merge in xyCoordM
        xyCoordM = [xyCoordM; [i timeMerge ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-1)+1:2*(timeMerge-1)+2) ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-2)+1:2*(timeMerge-2)+2) ...
            xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-2)+1:2*(timeMerge-2)+2)]];

        %add coordinates after merge to merging track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-1)+1:2*(timeMerge-1)+2) = ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-1)+1:2*(timeMerge-1)+2);

    end

    %get split events
    splitIndx = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)));

    %go over all splits
    for iSplit = splitIndx'

        %get split time
        timeSplit = seqOfEvents(iSplit,1);

        %store coordinates belonging to this split in xyCoordS
        xyCoordS = [xyCoordS; [i timeSplit ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2) ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-1)+1:2*(timeSplit-1)+2) ...
            xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-1)+1:2*(timeSplit-1)+2)]];

        %add coordinates before split to splitting track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-2)+1:2*(timeSplit-2)+2) = ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2);

    end

    %store coordinates in xyCoordAll
    xyCoordAll = [xyCoordAll; xyCoordTmp];

end

%replace NaNs by zero
xyCoordAll(isnan(xyCoordAll)) = 0;
