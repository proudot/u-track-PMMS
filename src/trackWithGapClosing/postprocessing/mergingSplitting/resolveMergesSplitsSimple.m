function tracksOut = resolveMergesSplitsSimple(tracksIn)
%RESOLVEMERGESSPLITSSIMPLE generates individual trajectories out of simple compound tracks 
%
%SYNPOSIS tracksOut = resolveMergesSplitsSimple(tracksIn)
%
%INPUT  tracksIn : Output of trackCloseGapsKalman.
%
%OUTPUT tracksOut: Similar to tracksIn, but with compound tracks resolved
%                  into individual tracks. In addition to tracks that
%                  exhibit no merging or splitting, only the following
%                  compound tracks are resolved: 
%                  (1) single merge
%                  (2) single split
%                  (3) split followed by merge
%                  (4) merge followed by split
%                  All other compound tracks are discarded.
%
%Khuloud Jaqaman, February 2009

%get number of tracks
numTracks = length(tracksIn);

%get number of segments in each track
numSegments = ones(numTracks,1);
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksIn(iTrack).tracksCoordAmpCG,1);
end

%keep only tracks with at most 3 segments
trackIndx123 = find(numSegments<=3);
tracksIn = tracksIn(trackIndx123);
numSegments = numSegments(trackIndx123);
numTracks = length(trackIndx123);

%get the indices of tracks with 1, 2 or 3 segments
trackIndx1 = find(numSegments==1);
numIndx1 = length(trackIndx1);
trackIndx2 = find(numSegments==2);
numIndx2 = length(trackIndx2);
trackIndx3 = find(numSegments==3);
numIndx3 = length(trackIndx3);

%initialize output structure and put tracks with only 1 segment in the
%output structure
tracksOut = repmat(struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],...
    'seqOfEvents',[]),2*numTracks,1);
if numIndx1 > 0
    tracksOut(1:numIndx1) = tracksIn(trackIndx1);
end

%go over tracks with 2 segments ...
%here there can be 3 scenarios: 
%(1) a single merge
%(2) a single split
%(3) a split followed by a merge
for iTrack = 1 : numIndx2

    %get current track index
    indxCurrent = trackIndx2(iTrack);

    %get current track's sequence of events
    seqOfEvents = tracksIn(indxCurrent).seqOfEvents;

    %resolve merges and splits based on the sequence of events
    if isnan(seqOfEvents(2,4)) %scenario (1) = single merge
        
        %get merge time relative to track start time
        mergeTime = seqOfEvents(3,1) - seqOfEvents(1,1) + 1;
        
        %copy track coordinates and amplitudes
        tracksCoordAmpCG = tracksIn(indxCurrent).tracksCoordAmpCG;
        
        %for the period after merging, give the same coordinates to both tracks
        newValues = max(tracksCoordAmpCG(1,(mergeTime-1)*8+1:end),tracksCoordAmpCG(2,(mergeTime-1)*8+1:end));
        tracksCoordAmpCG(:,(mergeTime-1)*8+1:end) = repmat(newValues,2,1);
        
        %find the ratio of amplitudes before merging, and distribute the
        %amplitude after merging accordingly
        ampTracks = tracksCoordAmpCG(:,4:8:(mergeTime-1)*8);
        ampTracks = nanmean(ampTracks,2);
        ampRatio1 = ampTracks(1)/(sum(ampTracks));
        tracksCoordAmpCG(1,(mergeTime-1)*8+4:8:end) = tracksCoordAmpCG(1,...
            (mergeTime-1)*8+4:8:end) * ampRatio1;
        tracksCoordAmpCG(2,(mergeTime-1)*8+4:8:end) = tracksCoordAmpCG(2,...
            (mergeTime-1)*8+4:8:end) * (1-ampRatio1);
        
        %write the new sequence of events
        seqOfEvents1 = [seqOfEvents(1,:); seqOfEvents(4,:)];
        seqOfEvents1(:,3) = 1;
        seqOfEvents2 = [seqOfEvents(2,:); seqOfEvents(4,:)];
        seqOfEvents2(:,3) = 1;
        secondTrackAppearance = seqOfEvents2(1,1) - seqOfEvents1(1,1) + 1;
        
        %separate the two track segments
        tracksCoordAmpCG1 = tracksCoordAmpCG(1,:);
        tracksCoordAmpCG2 = tracksCoordAmpCG(2,(secondTrackAppearance-1)*8+1:end);
        
        %store in output structure
        tracksOut(numIndx1+2*iTrack-1).tracksCoordAmpCG = tracksCoordAmpCG1;
        tracksOut(numIndx1+2*iTrack-1).seqOfEvents = seqOfEvents1;
        tracksOut(numIndx1+2*iTrack).tracksCoordAmpCG = tracksCoordAmpCG2;
        tracksOut(numIndx1+2*iTrack).seqOfEvents = seqOfEvents2;

    else %if second event is a split

        if isnan(seqOfEvents(3,4)) %scenario (2) = single split

            %get split time relative to track start time
            splitTime = seqOfEvents(2,1) - seqOfEvents(1,1) + 1;

            %copy track coordinates and amplitudes
            tracksCoordAmpCG = tracksIn(indxCurrent).tracksCoordAmpCG;

            %for the period before splitting, give the same coordinates to both tracks
            newValues = max(tracksCoordAmpCG(1,1:(splitTime-1)*8),tracksCoordAmpCG(2,1:(splitTime-1)*8));
            tracksCoordAmpCG(:,1:(splitTime-1)*8) = repmat(newValues,2,1);

            %find the ratio of amplitudes after splitting, and distribute the
            %amplitude before splitting accordingly
            ampTracks = tracksCoordAmpCG(:,(splitTime-1)*8+4:8:end);
            ampTracks = nanmean(ampTracks,2);
            ampRatio1 = ampTracks(1)/(sum(ampTracks));
            tracksCoordAmpCG(1,4:8:(splitTime-1)*8) = tracksCoordAmpCG(1,...
                4:8:(splitTime-1)*8) * ampRatio1;
            tracksCoordAmpCG(2,4:8:(splitTime-1)*8) = tracksCoordAmpCG(2,...
                4:8:(splitTime-1)*8) * (1-ampRatio1);

            %write the new sequence of events
            firstTrackEndRow = find(seqOfEvents(3:4,3)==1) + 2;
            secondTrackEndRow = setdiff((3:4),firstTrackEndRow);
            seqOfEvents1 = [seqOfEvents(1,:); seqOfEvents(firstTrackEndRow,:)];
            seqOfEvents1(:,3) = 1;
            seqOfEvents2 = [seqOfEvents(1,:); seqOfEvents(secondTrackEndRow,:)];
            seqOfEvents2(:,3) = 1;

            %separate the two track segments
            tracksCoordAmpCG1 = tracksCoordAmpCG(1,1:(seqOfEvents1(2,1)-seqOfEvents1(1,1)+1)*8);
            tracksCoordAmpCG2 = tracksCoordAmpCG(2,1:(seqOfEvents2(2,1)-seqOfEvents2(1,1)+1)*8);

            %store in output structure
            tracksOut(numIndx1+2*iTrack-1).tracksCoordAmpCG = tracksCoordAmpCG1;
            tracksOut(numIndx1+2*iTrack-1).seqOfEvents = seqOfEvents1;
            tracksOut(numIndx1+2*iTrack).tracksCoordAmpCG = tracksCoordAmpCG2;
            tracksOut(numIndx1+2*iTrack).seqOfEvents = seqOfEvents2;
            
        else %scenario (3) = split followed by merge

            %get split time relative to track start time
            splitTime = seqOfEvents(2,1) - seqOfEvents(1,1) + 1;

            %get merge time relative to track start time
            mergeTime = seqOfEvents(3,1) - seqOfEvents(1,1) + 1;

            %copy track coordinates and amplitudes
            tracksCoordAmpCG = tracksIn(indxCurrent).tracksCoordAmpCG;

            %for the period before splitting, give the same coordinates to both tracks
            newValues = max(tracksCoordAmpCG(1,1:(splitTime-1)*8),tracksCoordAmpCG(2,1:(splitTime-1)*8));
            tracksCoordAmpCG(:,1:(splitTime-1)*8) = repmat(newValues,2,1);

            %for the period after merging, give the same coordinates to both tracks
            newValues = max(tracksCoordAmpCG(1,(mergeTime-1)*8+1:end),tracksCoordAmpCG(2,(mergeTime-1)*8+1:end));
            tracksCoordAmpCG(:,(mergeTime-1)*8+1:end) = repmat(newValues,2,1);

            %find the ratio of amplitudes after splitting and before
            %merging, and distribute the amplitude before splitting and
            %after merging accordingly
            ampTracks = tracksCoordAmpCG(:,(splitTime-1)*8+4:8:(mergeTime-1)*8);
            ampTracks = nanmean(ampTracks,2);
            ampRatio1 = ampTracks(1)/(sum(ampTracks));
            tracksCoordAmpCG(1,4:8:(splitTime-1)*8) = tracksCoordAmpCG(1,...
                4:8:(splitTime-1)*8) * ampRatio1;
            tracksCoordAmpCG(2,4:8:(splitTime-1)*8) = tracksCoordAmpCG(2,...
                4:8:(splitTime-1)*8) * (1-ampRatio1);
            tracksCoordAmpCG(1,(mergeTime-1)*8+4:8:end) = tracksCoordAmpCG(1,...
                (mergeTime-1)*8+4:8:end) * ampRatio1;
            tracksCoordAmpCG(2,(mergeTime-1)*8+4:8:end) = tracksCoordAmpCG(2,...
                (mergeTime-1)*8+4:8:end) * (1-ampRatio1);

            %write the new sequence of events
            seqOfEvents1 = [seqOfEvents(1,:); seqOfEvents(4,:)];
            seqOfEvents1(:,3) = 1;
            seqOfEvents2 = seqOfEvents1;

            %separate the two track segments
            tracksCoordAmpCG1 = tracksCoordAmpCG(1,:);
            tracksCoordAmpCG2 = tracksCoordAmpCG(2,:);

            %store in output structure
            tracksOut(numIndx1+2*iTrack-1).tracksCoordAmpCG = tracksCoordAmpCG1;
            tracksOut(numIndx1+2*iTrack-1).seqOfEvents = seqOfEvents1;
            tracksOut(numIndx1+2*iTrack).tracksCoordAmpCG = tracksCoordAmpCG2;
            tracksOut(numIndx1+2*iTrack).seqOfEvents = seqOfEvents2;

        end %(if isnan(seqOfEvents(3,4)) ... else ...)

    end %(if isnan(seqOfEvents(2,4)) ... else ...)

end %(for iTrack = 1 : numIndx2)

%go over tracks with 3 segments
%here only accept one scenario, where 2 tracks merge then split
%the assignment of which track before merging goes to which track after
%splitting is based on intensity
iTrackUse = 0;
for iTrack = 1 : numIndx3

    %get current track index
    indxCurrent = trackIndx3(iTrack);

    %get current track's sequence of events
    seqOfEvents = tracksIn(indxCurrent).seqOfEvents;

    %resolve this compound track only if it falls in the following scenario
    if (seqOfEvents(2,2)==1 && isnan(seqOfEvents(2,4))) && ... %second event is an appearance
            (seqOfEvents(3,2)==2 && ~isnan(seqOfEvents(3,4))) && ... %third event is a merge
            (seqOfEvents(4,2)==1 && ~isnan(seqOfEvents(4,4))) && ... %fourth event is a split
            (seqOfEvents(5,2)==2 &&  isnan(seqOfEvents(5,4))) %fifth event is a disappearance
        
        %update the index iTrackUse
        iTrackUse = iTrackUse + 1;

        %get merge time relative to track start time
        mergeTime = seqOfEvents(3,1) - seqOfEvents(1,1) + 1;

        %get split time relative to track start time
        splitTime = seqOfEvents(4,1) - seqOfEvents(1,1) + 1;

        %copy track coordinates and amplitudes
        tracksCoordAmpCG = tracksIn(indxCurrent).tracksCoordAmpCG;

        %for the period after merging and before splitting, give the same
        %coordinates to all tracks
        newValues = max(tracksCoordAmpCG(1,(mergeTime-1)*8+1:(splitTime-1)*8),...
            tracksCoordAmpCG(2,(mergeTime-1)*8+1:(splitTime-1)*8));
        newValues = max(newValues,tracksCoordAmpCG(3,(mergeTime-1)*8+1:(splitTime-1)*8));
        tracksCoordAmpCG(:,(mergeTime-1)*8+1:(splitTime-1)*8) = repmat(newValues,3,1);

        %get average amplitudes before merging and indicate which segment
        %has smaller amplitude
        ampTracksBefore = tracksCoordAmpCG(1:2,4:8:(mergeTime-1)*8);
        ampTracksBefore = nanmean(ampTracksBefore,2);
        smallerAmpBeforeIndx = find(ampTracksBefore==min(ampTracksBefore));
        smallerAmpBeforeIndx = smallerAmpBeforeIndx(1);
        largerAmpBeforeIndx = setdiff((1:2),smallerAmpBeforeIndx);
        
        %get average amplitudes after splitting and indicate which segment
        %has smaller amplitude
        ampTracksAfter = tracksCoordAmpCG(:,(splitTime-1)*8+4:8:end);
        ampTracksAfter = nanmean(ampTracksAfter,2);
        smallerAmpAfterIndx = find(ampTracksAfter==min(ampTracksAfter));
        smallerAmpAfterIndx = smallerAmpAfterIndx(1);
        largerAmpAfterIndx = find(ampTracksAfter==max(ampTracksAfter));
        largerAmpAfterIndx = largerAmpAfterIndx(end);
        
        %extract the track data after the split, while distinguishing
        %between smaller and larger amplitude
        smallerAmpTrackData = tracksCoordAmpCG(...
            smallerAmpAfterIndx,(splitTime-1)*8+1:end);
        largerAmpTrackData = tracksCoordAmpCG(...
            largerAmpAfterIndx,(splitTime-1)*8+1:end);
        
        %put the smaller amplitude segments together, and the larger
        %amplitude segments together
        tracksCoordAmpCG(smallerAmpBeforeIndx,(splitTime-1)*8+1:end) = ...
            smallerAmpTrackData;
        tracksCoordAmpCG(largerAmpBeforeIndx,(splitTime-1)*8+1:end) = ...
            largerAmpTrackData;

        %get the amplitude ratio before merging and after splitting, and
        %distribute the amplitude after merging and before splitting
        %accordingly
        ampTracks = tracksCoordAmpCG(1:2,[4:8:(mergeTime-1)*8 (splitTime-1)*8+4:8:end]);
        ampTracks = nanmean(ampTracks,2);
        ampRatio1 = ampTracks(1)/(sum(ampTracks));
        tracksCoordAmpCG(1,(mergeTime-1)*8+4:8:(splitTime-1)*8) = tracksCoordAmpCG(1,...
            (mergeTime-1)*8+4:8:(splitTime-1)*8) * ampRatio1;
        tracksCoordAmpCG(2,(mergeTime-1)*8+4:8:(splitTime-1)*8) = tracksCoordAmpCG(2,...
            (mergeTime-1)*8+4:8:(splitTime-1)*8) * (1-ampRatio1);
        
        %find the track end times
        trackSEL = getTrackSEL(tracksCoordAmpCG);
        track1End = trackSEL(1,2);
        track2End = trackSEL(2,2);

        %write the new sequence of events
        seqOfEvents1 = [seqOfEvents(1,:); [track1End 2 1 NaN]];
        seqOfEvents1(:,3) = 1;
        seqOfEvents2 = [seqOfEvents(2,:); [track2End 2 1 NaN]];
        seqOfEvents2(:,3) = 1;
        secondTrackAppearance = seqOfEvents2(1,1) - seqOfEvents1(1,1) + 1;

        %separate the two track segments
        tracksCoordAmpCG1 = tracksCoordAmpCG(1,1:(seqOfEvents1(2,1)-seqOfEvents1(1,1)+1)*8);
        tracksCoordAmpCG2 = tracksCoordAmpCG(2,(secondTrackAppearance-1)*8+1:(seqOfEvents2(2,1)-seqOfEvents1(1,1)+1)*8);

        %store in output structure
        tracksOut(numIndx1+2*numIndx2+2*iTrackUse-1).tracksCoordAmpCG = tracksCoordAmpCG1;
        tracksOut(numIndx1+2*numIndx2+2*iTrackUse-1).seqOfEvents = seqOfEvents1;
        tracksOut(numIndx1+2*numIndx2+2*iTrackUse).tracksCoordAmpCG = tracksCoordAmpCG2;
        tracksOut(numIndx1+2*numIndx2+2*iTrackUse).seqOfEvents = seqOfEvents2;

    end %(if (seqOfEvents(2,2)==1 && ...)

end %(for iTrack = 1 : numIndx3)

tracksOut = tracksOut(1:numIndx1+2*numIndx2+2*iTrackUse);

% ~~~ the end ~~~