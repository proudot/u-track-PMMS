function [linkStats,gapStats] = scoreLinksGaps(tracksFinal,simMPM)
%SCORELINKSGAPS compares links and closed gaps obtained via tracking to the ground truth
%
%SYNOPSIS [linkStats,gapStats] = scoreLinksGaps(tracksFinal,simMPM)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix).
%       simMPM      : Simulated tracks, as obtained from
%                     simulateMimickCD36 or simulateMPMlftbeh.
%
%OUTPUT linkStats   : (number of frames - 1) - by - 5 array. Row i 
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
%       
%Khuloud Jaqaman, September 2007

%% convert tracksFinal from structure to matrix

%get number of frames
numFrames = size(simMPM,2)/3;

%convert structure to matrix if necessary
if isstruct(tracksFinal) %if tracks are input in structure format

    %get number of tracks and frames
    numTracks = length(tracksFinal);
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
    end

    %locate the row of the first track of each compound track in the
    %big matrix of all tracks (to be constructed in the next step)
    trackStartRow = ones(numTracks,1);
    for iTrack = 2 : numTracks
        trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
    end

    %put all tracks together in a matrix
    trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numFrames);
    for i = 1 : numTracks
        startTime = tracksFinal(i).seqOfEvents(1,1);
        endTime   = tracksFinal(i).seqOfEvents(end,1);
        trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
            numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
            tracksFinal(i).tracksCoordAmpCG;
    end
    
else %if tracks are input in matrix format
    
    %copy input into new variable
    trackedFeatureInfo = tracksFinal;
    
end
clear tracksFinal

%retain only the x and y coordinates from trackedFeatureInfo and simMPM
xCoord1 = trackedFeatureInfo(:,1:8:end);
xCoord1(isnan(xCoord1)) = 0;
yCoord1 = trackedFeatureInfo(:,2:8:end);
yCoord1(isnan(yCoord1)) = 0;
xCoord0 = simMPM(:,1:3:end);
yCoord0 = simMPM(:,2:3:end);

%% establish correspondence between features in tracking results and in ground truth

%pre-allocate memory
tracksCorrespond = zeros(size(xCoord1));
tracksGT = tracksCorrespond;

%establish feature correspondence between tracking results and ground truth
for iFrame = 1 : numFrames
    
    %get ground truth feature coordinates in current frame
    xyCoord0 = [yCoord0(:,iFrame) xCoord0(:,iFrame)];
    xyCoord0 = xyCoord0(xyCoord0(:,1)~=0,:);
    numFeat0 = size(xyCoord0,1);
    
    %get tracking results feature coordinates in current frame
    xyCoord1 = [xCoord1(:,iFrame) yCoord1(:,iFrame)];
    xyCoord1 = xyCoord1(xyCoord1(:,1)~=0,:);
    numFeat1 = size(xyCoord1,1);
    
    %calculate the distance matrix between ground truth features and
    %tracking results features
    distanceMat = createDistanceMatrix(xyCoord1,xyCoord0);
    distanceMat(distanceMat>2) = -1; %assume that detection won't be wrong by more than 2 pixels
    
    %ironically, use LAP to get correspondence
    link10 = lap(distanceMat,-1,0,1);
    
    %store correspondence
    link10 = link10(1:numFeat1);
    link10(link10>numFeat0) = -1; %false positives that have no counterpart in ground truth
    tracksCorrespond(xCoord1(:,iFrame)~=0,iFrame) = link10;
    
    %store indices of features in ground truth
    tracksGT(xCoord0(:,iFrame)~=0,iFrame) = (1:numFeat0)';
        
end

%% compare links from tracking code to ground truth links

%initialize total number of links and number of correct links
linkStats = NaN(numFrames-1,5);

%go over all frames ...
for iFrame = 1 : numFrames - 1
    
    %% number of links in ground truth - all correct by definition ...
    
    %get indices in ground truth
    tracksGT12 = tracksGT(:,iFrame:iFrame+1);
    tracksGT12 = tracksGT12(tracksGT12(:,1)~=0,:);
    
    %get number of all links in ground truth
    numLinksGT12 = length(find( tracksGT(:,1) ~= 0 & tracksGT(:,2) ~= 0 ));
    
    %% number of links from tracking code - some are correct, some are
    %% between real features but are wrong, others involve a false
    %% detection positive and thus are wrong ...
    
    %get the indices of features in ground truth corresponding to features
    %in tracking code results in iFrame and iFrame+1
    tracksCorrespond12 = tracksCorrespond(:,iFrame:iFrame+1);
    
    %remove zeros from 1st frame
    tracksCorrespond12 = tracksCorrespond12(tracksCorrespond12(:,1)~=0,:);
    
    %get number of all links from the tracking code
    numLinksTrackCode12 = length(find( tracksCorrespond12(:,1) ~= 0 & ...
        tracksCorrespond12(:,2) ~= 0 ));
    
    %% number of correct links from tracking code ...
    
    %find features that are in common in the first frame between ground
    %truth and tracking code results (i.e., which features in the ground 
    %truth survived the detection step)
    [dummy,ia,ib] = intersect(tracksGT12(:,1),tracksCorrespond12(:,1));
    
    %compare the features they get linked to in both tracking results and
    %ground truth - those that are equal are correct links
    numLinksCorrect12 = length(find( tracksGT12(ia,2)==tracksCorrespond12(ib,2) ...
        & tracksGT12(ia,2)~=0 ));
    
    %% number of wrong links from tracking code that are between real features ...
    
    %get number of all links from tracking code between
    %real features (i.e. ignoring false positives from the detector)
    numLinksReal12 = length(find( tracksCorrespond12(ib,2)~=0 & ...
        tracksCorrespond12(ib,2)~=-1 ));
    
    %hence calculate number of wrong links between real features
    numLinksWrong12 = numLinksReal12 - numLinksCorrect12;
    
    %% number of links from tracking code involving false detection
    %% positives ...
        
    %get number of all links in tracking results that involve a false
    %positive (which are obviously wrong)
    numLinksFalse12 = length(find( (tracksCorrespond12(:,1) == -1 & ...
        tracksCorrespond12(:,2) ~= 0) | (tracksCorrespond12(:,2) == -1) ));
    
    %% store numbers for output ...
    
    linkStats(iFrame,:) = [numLinksGT12 numLinksTrackCode12 numLinksCorrect12 ...
        numLinksWrong12 numLinksFalse12];
    
end

%% evaluate gaps

%find gaps in tracks
gapInfo = findTrackGaps(trackedFeatureInfo);
numGaps = size(gapInfo,1);

%allocate memory for gapStats
gapStats = NaN(numGaps,2);

%go over all gaps ...
for iGap = 1 : numGaps
    
    %get some gap information
    iTrack = gapInfo(iGap,1);
    frameBef = gapInfo(iGap,3) - 1;
    gapLength = gapInfo(iGap,4);
    frameAft = frameBef + gapLength + 1;
    
    %find the indices of ground truth features corresponding to the
    %features from linking results
    featBef1 = tracksCorrespond(iTrack,frameBef);
    featAft1 = tracksCorrespond(iTrack,frameAft);
    
    if featBef1 ~= -1 && featAft1 ~= -1
        
        %find featBef1 in frameBef in the ground truth tracks and look up which
        %feature it links to in the ground truth in frameAft
        featAft0 = tracksGT(tracksGT(:,frameBef)==featBef1,frameAft);

        %store wether gap is correctly closed (0) or not (1)
        gapStats(iGap,:) = [gapLength ~(featAft0==featAft1)];

    else

        %indicate that gap involves detection artifacts
        gapStats(iGap,:) = [gapLength 2];
        
    end
        
end



%% ~~~ the end ~~~ %%

