function linkStatsLag = scoreLinksWithLag(tracksFinal,tracksSim)
%SCORELINKSWITHLAG compares links obtained via tracking to the ground truth, looking at different lags
%
%SYNOPSIS linkStatsLag = scoreLinksWithLag(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : Simulated tracks as obtained from
%                     simulateMimickCD36_MS (structure format).
%
%OUTPUT linkStatsLag: (number of frames - 1) - by - 5 array. Row i 
%                     corresponds to the links from frame t to frame t+i.
%                     The columns show the number of:
%                     (1) ground truth links;
%                     (2) links resulting from tracking;
%                     (3) tracking links that are correct;
%                     (4) tracking links that are between real features but
%                         are wrong;
%       
%Khuloud Jaqaman, October 2007

%% process input variables

%extract track information out of tracksSim
xyCoordAll0 = trackInfoFromStruct(tracksSim);

if isstruct(tracksFinal)

    %extract track information out of tracksFinal
    [xyCoordAll1,dummy,dummy,trackStartRow1] = ...
        trackInfoFromStruct(tracksFinal);

else

    %directly extract track information out of tracksFinal
    xyCoordAll1 = zeros(size(tracksFinal,1),size(tracksFinal,2)/4);
    xyCoordAll1(:,1:2:end) = tracksFinal(:,1:8:end);
    xyCoordAll1(:,2:2:end) = tracksFinal(:,2:8:end);
    xyCoordAll1(isnan(xyCoordAll1)) = 0;
    trackStartRow1 = (1:size(xyCoordAll1,1))';
    
end

%get number of frames in ground truth and in tracking results
numFrames0 = size(xyCoordAll0,2)/2;
numFrames1 = size(xyCoordAll1,2)/2;
if numFrames0 ~= numFrames1
    numFrames0 = min(numFrames0,numFrames1);
    disp('ATTENTION: different number of frames in ground truth and simulation results!')
end

%% evaluate gaps and store gap status in matrix of coordinates

%find gaps in tracks from tracking code
gapInfo = findTrackGaps(tracksFinal);
numGaps = size(gapInfo,1);

%go over all gaps ...
for iGap = 1 : numGaps
    
    %get some gap information
    iTrack = gapInfo(iGap,1);
    iSegment = gapInfo(iGap,2);
    frameBef = gapInfo(iGap,3) - 1;
    gapLength = gapInfo(iGap,4);
    frameAft = frameBef + gapLength + 1;
    
    %find the coordinates from the tracking results before and after the gap
    xyCoordBef1 = xyCoordAll1(trackStartRow1(iTrack)+iSegment-1,(frameBef-1)*2+1:(frameBef-1)*2+2);
    xyCoordAft1 = xyCoordAll1(trackStartRow1(iTrack)+iSegment-1,(frameAft-1)*2+1:(frameAft-1)*2+2);
    
    %determine whether either feature is a false positive
    realFeatBef = ismember(xyCoordBef1,xyCoordAll0(:,(frameBef-1)*2+1:(frameBef-1)*2+2),'rows');
    realFeatAft = ismember(xyCoordAft1,xyCoordAll0(:,(frameAft-1)*2+1:(frameAft-1)*2+2),'rows');
    
    if realFeatBef && realFeatAft
        
        %find featBef1 in frameBef in the ground truth tracks and look up which
        %feature it links to in the ground truth in frameAft
        indx = find( (xyCoordBef1(1) == xyCoordAll0(:,(frameBef-1)*2+1)) & ...
            (xyCoordBef1(2) == xyCoordAll0(:,(frameBef-1)*2+2)),1,'first' );
        
        xyCoordAft0 = xyCoordAll0(indx,(frameAft-1)*2+1:(frameAft-1)*2+2);

        %store wether gap is correctly closed (0) or not (1)
        gapInfo(iGap,7) = ~min((xyCoordAft1 == xyCoordAft0));

    else

        %indicate that gap involves detection artifacts
        gapInfo(iGap,7) = 2;
        
    end
    
    %store status of gap in xyCoordAll1
    xyCoordAll1(trackStartRow1(iTrack)+iSegment-1,frameBef*2+1:(frameAft-1)*2) = -1 - gapInfo(iGap,7);
        
end

%% score links

%initialize output variable
linkStatsLag = NaN(numFrames0-1,4);

%go over all lags ...
for lag = 1 : numFrames0 - 1

    %initialize linkStats
    linkStats = zeros(1,4);
    
    %go over all relevant frames ...
    for iFrame = 1 : numFrames0 - lag

        %% number of links in ground truth - all correct by definition ...

        %get coordinates in ground truth
        tracksGT12 = [xyCoordAll0(:,iFrame*2-1:iFrame*2) ...
            xyCoordAll0(:,(iFrame+lag)*2-1:(iFrame+lag)*2)];

        %get number of links in ground truth
        numLinksGT12 = length(find( tracksGT12(:,1) ~= 0 & tracksGT12(:,3) ~= 0 ));

        %% number of links from tracking code - some are correct, some are
        %% between real features but are wrong, others involve a false
        %% detection positive and thus are wrong ...

        %get coordinates in tracking results
        tracksCode12 = [xyCoordAll1(:,iFrame*2-1:iFrame*2) ...
            xyCoordAll1(:,(iFrame+lag)*2-1:(iFrame+lag)*2)];

        %get number of links in tracking results ("links" to "gaps" are included)
        numLinksTrackCode12 = length(find( tracksCode12(:,1)~=0 & tracksCode12(:,3)~=0 ));

        %% number of correct links from tracking code (links between features and links to gaps) ...

        %get links between features in tracking results which are found in the ground truth
        commonLinks = ismember(tracksCode12,tracksGT12,'rows');

        %remove zeros to get number of correct links between features
        numLinksCorrect12 = length(find( commonLinks(:,1)~=0 & ...
            tracksCode12(:,1)~=0 & tracksCode12(:,3)~=0 ));
        
        %add links that include a correctly closed gap
        numLinksCorrect12 = numLinksCorrect12 + ...
            length(find( ismember(tracksCode12(:,1:2),tracksGT12(:,1:2),'rows') & ...
            tracksCode12(:,1)~=0 & tracksCode12(:,3)==-1 )) + ...
            length(find( ismember(tracksCode12(:,3:4),tracksGT12(:,3:4),'rows') & ...
            tracksCode12(:,3)~=0 & tracksCode12(:,1)==-1 )) + ...
            length(find( tracksCode12(:,1)==-1 & tracksCode12(:,3)==-1 ));

        %% number of wrong links from tracking code that are between real features ...

        %find the "real" features in the 2 frames
        realFeat1 = ismember(tracksCode12(:,1:2),tracksGT12(tracksGT12(:,1)~=0,1:2),'rows');
        realFeat2 = ismember(tracksCode12(:,3:4),tracksGT12(tracksGT12(:,3)~=0,3:4),'rows');

        %get number of all links from tracking code between real features
        numLinksReal12 = length(find( realFeat1~=0 & realFeat2~=0 ));
        
        %add links that involve gaps (whether correctly or wrongly closed, but not involving false positives)
        numLinksReal12 = numLinksReal12 + ...
            length(find( realFeat1~=0 & tracksCode12(:,3)<0 & tracksCode12(:,3)~=-3 )) + ...
            length(find( realFeat2~=0 & tracksCode12(:,1)<0 & tracksCode12(:,1)~=-3 )) + ...
            length(find( tracksCode12(:,1)<0 & tracksCode12(:,1)~=-3 & ...
            tracksCode12(:,3)<0 & tracksCode12(:,3)~=-3 ));

        %hence calculate number of wrong links between real features
        numLinksWrong12 = numLinksReal12 - numLinksCorrect12;            
        
        %% store numbers for output ...

        linkStats = linkStats + [numLinksGT12 numLinksTrackCode12 ...
            numLinksCorrect12 numLinksWrong12];

    end %(for iFrame = 1 : numFrames0 - lag)
    
    %store linking statistics for this lag
    linkStatsLag(lag,:) = linkStats;
    
end %(for lag = 1 : numFrames0 - 1)

%% Subfunction 1

function [xyCoordAll,xyCoordM,xyCoordS,trackStartRow] = trackInfoFromStruct(tracks)

%get number of tracks
numTracks = length(tracks);

%get number of frames
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%find gaps in tracks
gapInfo = findTrackGaps(tracks);

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

    %fill in -1 for gaps
    indx = (find(gapInfo(:,1)==i))';
    for iGap = indx
        xyCoordTmp(gapInfo(iGap,2),gapInfo(iGap,3)*2-1:(gapInfo(iGap,3)+gapInfo(iGap,4)-1)*2) = -1;
    end
        
    %determine row where this compound track starts in xyCoordAll
    trackStartRow(i) = size(xyCoordAll,1) + 1;

    %get merge events
    mergeIndx = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));

    %go over all merges
    for iMerge = mergeIndx'

        %get merge time
        timeMerge = seqOfEvents(iMerge,1);

        %store coordinates belonging to this merge in xyCoordM
        xyCoordM = [xyCoordM; [timeMerge ...
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
        xyCoordS = [xyCoordS; [timeSplit ...
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
