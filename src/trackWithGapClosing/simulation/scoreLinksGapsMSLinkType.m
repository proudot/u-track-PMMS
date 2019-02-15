function [linkStats,gapStats,mergeSplitStats] = scoreLinksGapsMSLinkType(tracksFinal,tracksSim)
%SCORELINKSGAPSMS compares links, closed gaps and merges/splits obtained via tracking to the ground truth
%
%SYNOPSIS [linkStats,gapStats,mergeSplitStats] = scoreLinksGapsMS(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : Simulated tracks as obtained from
%                     simulateMimickCD36_MS (structure format).
%
%OUTPUT linkStats   : (Number of frames - 1) - by - 5 array. Row i 
%                     corresponds to the links from frame i to frame i+1.
%                     The columns show the number of:
%                     (1) ground truth links;
%                     (2) links resulting from tracking;
%                     (3) tracking links that are correct;
%                     (4) tracking links that are between real features but
%                         are wrong;
%                     (5) tracking links that involve a detection artifact
%                         (i.e. a detetion false positive).
%                     The sum of columns 3 through 5 = column 2.
%       gapStats    : Two-column array with number of rows = number of
%                     gaps. First column shows gap length. Second column is
%                     0 - if gap is correctly closed;
%                     1 - if features just before and just after the gap
%                     are real but the connection is wrong;
%                     2 - if one or both features just before and just
%                     after the gap are detection artifacts.
%       mergeSplitStats: 2-by-3 matrix where first row is for merges, 2nd
%                        row is for splits.  1st column shows number of
%                        merges/splits in ground truth, 2nd column shows
%                        number of merges/splits in tracking results, and
%                        3rd column shows number of correct merges/splits
%                        in tracking results. Merge/split statistics
%                        don't explicitly account for detection false
%                        positives, i.e. they simply count toward wrong
%                        merges/splits, just like real features if wrongly
%                        assigned a merge/split.
%       
%Khuloud Jaqaman, October 2007

%% process input variables

%extract track information out of tracksSim
[xyCoordAll0,xyCoordM0,xyCoordS0,dummy,linkType] = trackInfoFromStruct(tracksSim);

if isstruct(tracksFinal)

    %extract track information out of tracksFinal
    [xyCoordAll1,xyCoordM1,xyCoordS1,trackStartRow1] = ...
        trackInfoFromStruct(tracksFinal);

    %look at merges and splits
    mergeSplit = 1;

else

    %directly extract track information out of tracksFinal
    xyCoordAll1 = zeros(size(tracksFinal,1),size(tracksFinal,2)/4);
    xyCoordAll1(:,1:2:end) = tracksFinal(:,1:8:end);
    xyCoordAll1(:,2:2:end) = tracksFinal(:,2:8:end);
    xyCoordAll1(isnan(xyCoordAll1)) = 0;
    trackStartRow1 = (1:size(xyCoordAll1,1))';
    
    %do not look at merges and splits
    mergeSplit = 0;

end

%get number of frames in ground truth and in tracking results
numFrames0 = size(xyCoordAll0,2)/2;
numFrames1 = size(xyCoordAll1,2)/2;
if numFrames0 ~= numFrames1
    numFrames0 = min(numFrames0,numFrames1);
    disp('ATTENTION: different number of frames in ground truth and simulation results!')
end

%% compare links (implicitly including merges and splits)

%initialize output variable (4 dimension for 4 type.
linkStats = NaN(numFrames0-1,5,4);

% for every track type (0 is for every track)
for type =0:1:3

    %go over all frames ...
    for iFrame = 1 : numFrames0 - 1
        
        %% number of links in ground truth - all correct by definition ...
        
        %get coordinates in ground truth (four coordinates at the same time)
        if(type==0)
            tracksGT12 = xyCoordAll0(:,(iFrame-1)*2+1:(iFrame+1)*2);
        else
            tracksGT12= zeros(size(xyCoordAll0,1),4);
            tracksGT12(linkType(:,iFrame)==type,:) = xyCoordAll0(linkType(:,iFrame)==type,(iFrame-1)*2+1:(iFrame+1)*2);
        end 
        
        %get number of links in ground truth
        numLinksGT12 = length(find( tracksGT12(:,1) ~= 0 & tracksGT12(:,3) ~= 0 ));
        
        %% number of links from tracking code - some are correct, some are
        %% between real features but are wrong, others involve a false
        %% detection positive and thus are wrong ...
        
        %get coordinates in tracking results
        tracksCode12 = xyCoordAll1(:,(iFrame-1)*2+1:(iFrame+1)*2);
        
        %get number of links in tracking results
        numLinksTrackCode12 = length(find( tracksCode12(:,1)~=0 & tracksCode12(:,3)~=0 ));
        
        %% number of correct links from tracking code ...
        
        %get common links between ground truth and tracking results
        commonLinks = ismember(tracksCode12,tracksGT12,'rows');
        
        %remove zeros to get number of correct links
        numLinksCorrect12 = length(find( commonLinks~=0 & tracksCode12(:,1)~=0 & tracksCode12(:,3)~=0 ));

        %% number of wrong links from tracking code that are between real features ...
        
        %find the "real" features in the 2 frames
        realFeat1 = ismember(tracksCode12(:,1:2),tracksGT12(tracksGT12(:,1)~=0,1:2),'rows');
        realFeat2 = ismember(tracksCode12(:,3:4),tracksGT12(tracksGT12(:,3)~=0,3:4),'rows');
        
        %get number of all links from tracking code between real features
        numLinksReal12 = length(find( realFeat1~=0 & realFeat2~=0 ));
        
        %hence calculate number of wrong links between real features
        numLinksWrong12 = numLinksReal12 - numLinksCorrect12;
        
        %% number of links from tracking code involving false detection
        %% positives ...
        
        %get number of all links in tracking results that involve a false
        %positive (which are obviously wrong)
        numLinksFalse12 = length(find( ...
            (realFeat1 == 0 & tracksCode12(:,1) ~= 0 & tracksCode12(:,3) ~= 0) | ...
            (realFeat2 == 0 & tracksCode12(:,3) ~= 0 & tracksCode12(:,1) ~= 0) ));
        
        %% store numbers for output ...
        
        linkStats(iFrame,:,type+1) = [numLinksGT12 numLinksTrackCode12 numLinksCorrect12 ...
                            numLinksWrong12 numLinksFalse12];
        
    end
end
% $$$ linkStatsBrown=linkStats(:,:,2);        
% $$$ linkStatsConf=linkStats(:,:,3);
% $$$ linkStatsDirect=linkStats(:,:,4);
% $$$ linkStats=linkStats(:,:,1);

%% evaluate gaps

%find gaps in tracks from tracking code
gapInfo = findTrackGaps(tracksFinal);
numGaps = size(gapInfo,1);

%allocate memory for gapStats
gapStats = NaN(numGaps,2);

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
            (xyCoordBef1(2) == xyCoordAll0(:,(frameBef-1)*2+2)));
        
        xyCoordAft0 = xyCoordAll0(indx,(frameAft-1)*2+1:(frameAft-1)*2+2); %#ok<FNDSB>

        %store wether gap is correctly closed (0) or not (1)
        gapStats(iGap,:) = [gapLength ~ismember(xyCoordAft1,xyCoordAft0,'rows')];

    else

        %indicate that gap involves detection artifacts
        gapStats(iGap,:) = [gapLength 2];
        
    end
        
end

%% check merges and splits

if mergeSplit
    
    %initialize number of correct merges
    numCorrectMerges = 0;
    
    %get number of merges in ground truth and tracking results
    numMerge0 = size(xyCoordM0,1);
    numMerge1 = size(xyCoordM1,1);
        
    %go over all merges in tracking results ...
    for iMerge = 1 : numMerge1
        
        %get time of merge
        timeMerge = xyCoordM1(iMerge,1);
        
        %get coordinates before and after merge
        xyCoordMAft1  = xyCoordM1(iMerge,2:3);
        xyCoordMBef11 = xyCoordM1(iMerge,4:5);
        xyCoordMBef12 = xyCoordM1(iMerge,6:7);
        
        %check if this merge exists in the ground truth
        [~,indxGT] = ismember([timeMerge xyCoordMAft1],xyCoordM0(:,1:3),'rows');
        
        %if merge exists in ground truth ...
        if indxGT > 0
            
            %get coordinates of merging features in ground truth
            xyCoordMBef01 = xyCoordM0(indxGT,4:5);
            xyCoordMBef02 = xyCoordM0(indxGT,6:7);
            
            %find how many tracking and ground truth coordinates are equivalent
            numEquiv = sum( ismember([xyCoordMBef11; xyCoordMBef12],...
                [xyCoordMBef01; xyCoordMBef02],'rows') );
            
            %update number of correct merges based on numEquiv
            switch numEquiv
                
                case 2 %if both coordinates have equivalents
                    
                    %add one to number of correct merges
                    numCorrectMerges = numCorrectMerges + 1;
                    
                case 1 %if only one coordinate has an equivalent
                    
                    %add 1/2 to number of correct merges
                    numCorrectMerges = numCorrectMerges + 0.5;
                    
                case 0 %if neither coordinate has an equivalent
                    
                    %don't add anything to number of correct merges
                    
            end %(switch numEquiv)
            
        end %(if ~isempty(indxGT))
        
    end %(for iMerge = 1 : numMerge1)
    
    %initialize number of correct splits
    numCorrectSplits = 0;
    
    %get number of splits in ground truth and tracking results
    numSplit0 = size(xyCoordS0,1);
    numSplit1 = size(xyCoordS1,1);
        
    %go over all splits in tracking results ...
    for iSplit = 1 : numSplit1
        
        %get time of split
        timeSplit = xyCoordS1(iSplit,1);
        
        %get coordinates before and after split
        xyCoordSBef1  = xyCoordS1(iSplit,2:3);
        xyCoordSAft11 = xyCoordS1(iSplit,4:5);
        xyCoordSAft12 = xyCoordS1(iSplit,6:7);
        
        %check if this split exists in the ground truth
        [~,indxGT] = ismember([timeSplit xyCoordSBef1],xyCoordS0(:,1:3),'rows');
        
        %if split exists in ground truth ...
        if indxGT > 0
            
            %get coordinates of splitting features in ground truth
            xyCoordSAft01 = xyCoordS0(indxGT,4:5);
            xyCoordSAft02 = xyCoordS0(indxGT,6:7);
            
            %find how many tracking and ground truth coordinates are equivalent
            numEquiv = sum( ismember([xyCoordSAft11; xyCoordSAft12],...
                [xyCoordSAft01; xyCoordSAft02],'rows') );
            
            %update number of correct splits based on numEquiv
            switch numEquiv
                
                case 2 %if both coordinates have equivalents
                    
                    %add one to number of correct splits
                    numCorrectSplits = numCorrectSplits + 1;
                    
                case 1 %if only one coordinate has an equivalent
                    
                    %add 1/2 to number of correct splits
                    numCorrectSplits = numCorrectSplits + 0.5;
                    
                case 0 %if neither coordinate has an equivalent
                    
                    %don't add anything to number of correct splits
                    
            end %(switch numEquiv)
            
        end %(if ~isempty(indxGT))
        
    end %(for iSplit = 1 : numSplit1)
    
    mergeSplitStats = [numMerge0 numMerge1 numCorrectMerges; ...
        numSplit0 numSplit1 numCorrectSplits];
    
else
    
    mergeSplitStats = [];
    
end %(if mergeSplit)

%% Subfunction 1

function [xyCoordAll,xyCoordM,xyCoordS,trackStartRow,linkType] = trackInfoFromStruct(tracks)

%get number of tracks
numTracks = length(tracks);

%get number of frames
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%initialize output variables
xyCoordAll = [];
linkType = [];
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

    linkTypeTmp = zeros(size(tracksCoordAmpCG,1),numFrames);
    if(isfield(tracks(i),'tracksSegType'))
        linkTypeTmp(:,startTime:endTime)=tracks(i).tracksSegType;
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
            xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-2)+1:2*(timeMerge-2)+2)]]; %#ok<AGROW>
        
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
            xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-1)+1:2*(timeSplit-1)+2)]]; %#ok<AGROW>
        
        %add coordinates before split to splitting track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-2)+1:2*(timeSplit-2)+2) = ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2);

    end

    %store coordinates in xyCoordAll
    xyCoordAll = [xyCoordAll; xyCoordTmp]; %#ok<AGROW>
    linkType = [linkType; linkTypeTmp]; %#ok<AGROW>

end

%replace NaNs by zero
xyCoordAll(isnan(xyCoordAll)) = 0;
