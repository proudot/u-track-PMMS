function [costMat,noLinkCost,nonlinkMarker,trackStartTime,trackEndTime,...
    indxMerge,numMerge,indxSplit,numSplit,errFlag] = costMatCloseGaps(...
    trackedFeatInfo,trackStartTime,indxStart,trackEndTime,indxEnd,...
    costMatParams,gapCloseParam)
%COSTMATCLOSEGAPS provides a cost matrix for closing gaps (including merging and splitting) based on tracked feature statistics
%
%SYNOPSIS [costMat,noLinkCost,nonlinkMarker,trackStartTime,trackEndTime,...
%    indxMerge,numMerge,indxSplit,numSplit,errFlag] = costMatCloseGaps(...
%    trackedFeatInfo,trackStartTime,indxStart,trackEndTime,indxEnd,...
%    costMatParams,gapCloseParam)
%
%INPUT  trackedFeatInfo: The positions and amplitudes of the tracked
%                        features from linkFeaturesTp2Tp. 
%                        Number of rows = number of tracks.
%                        Number of columns = 8*number of time points. 
%                        Each row consists of 
%                        [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                        in image coordinate system (coordinates in
%                        pixels). NaN is used to indicate time points 
%                        where the track does not exist.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       costMatParams  : Structure with the fields:
%             .trackStats  : Structure with the following fields:
%                   .dispSqR     : timeWindow x 1 vector of r in the gamma
%                                  distribution that describes the displacement
%                                  of a feature between frames.
%                   .dispSqTheta : timeWindow x 1 vector of theta in the 
%                                  gamma distribution that describes the 
%                                  displacement of a feature between frames.
%                   .ampDiffStd  : Standard deviations of the change in a feature's 
%                                  amplitude between two time points.
%             .cutoffProbD1: Cumulative probability of a square diplacement
%                            beyond which linking between an end and a 
%                            start is not allowed.
%             .cutoffProbA1: Cumulative probability of an amplitude change
%                            beyond which linking between an end and a 
%                            start is not allowed.
%             .cutoffProbD2: Cumulative probability of a square displacement
%                            beyond which merging and splitting are not
%                            allowed.
%             .cutoffProbA2: Cumulative probability of an amplitude change
%                            beyond which merging and splitting are not
%                            allowed.
%             .noLnkPrctl  : Percentile used to calculate the cost of
%                            linking a feature to nothing. Use -1 if you do
%                            not want to calculate this cost.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to it.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%
%OUTPUT costMat       : Cost matrix.
%       noLinkCost    : Cost of linking a feature to nothing, as derived
%                       from the distribution of costs.
%       nonlinkMarker : Value indicating that a link is not allowed.
%       trackStartTime: Starting time of tracks whose starts are considered
%                       for gap closing.
%       trackEndTime  : Ending time of tracks whose ends are considered for
%                       gap closing.
%       indxMerge     : Index of tracks that have possibly merged with
%                       tracks that end before the last time points.
%       numMerge      : Number of such tracks.
%       indxSplit     : Index of tracks from which tracks that begin after
%                       the first time point might have split.
%       numSplit      : Number of such tracks.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%REMARKS the costs are given by ...
%
%The cost for linking the end of track i at time point t to the 
%beginning of track j at time point t' (t' > t) is given by
%-log[p(dI)p(dispSq)], where dI = I(j;t') - I(i,t) is assumed to be 
%normally distributed with mean 0 and standard deviation ampDiffStd(t'-t) 
%(supplied by user), and dispSq = square of distance between end position 
%of i at t and beginning position of j at t' is assumed to be gamma
%distributed with parameters dispSqR(t'-t) and dispSqTheta(t'-t) (supplied by user).
%
%The cost for merging and splitting is given by -log[p(dI)p(dispSq)], where
%dispSq and p(dispSq) are as described above and, in the case of merging,
%dI = I(feature after merging) - sum(I(features before merging)), and, 
%in the case of splitting, 
%dI = sum(I(features after splitting)) - I(feature before splitting). 
%dI is thus normally distributed with mean zeros and standard deviation
%sqrt(2)*ampDiffStd(1).
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costMat = [];
noLinkCost = [];
nonlinkMarker = [];
indxMerge = [];
numMerge = [];
indxSplit = [];
numSplit = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatCloseGaps')
    disp('--costMatCloseGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

timeWindow = gapCloseParam.timeWindow;
mergeSplit = gapCloseParam.mergeSplit;
dispSqR = costMatParams.trackStats.dispSqR;
dispSqTheta = costMatParams.trackStats.dispSqTheta;
ampDiffStd = costMatParams.trackStats.ampDiffStd;
cutoffProbD1 = costMatParams.cutoffProbD1;
cutoffProbA1 = costMatParams.cutoffProbA1;
if mergeSplit
    cutoffProbD2 = costMatParams.cutoffProbD2;
    cutoffProbA2 = costMatParams.cutoffProbA2;
end
noLnkPrctl = costMatParams.noLnkPrctl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the maximum squared displacement that allows linking an end to a start
maxDispSq = expinv(cutoffProbD1,dispSqR,1./dispSqTheta);

%find the maximum squared amplitude change that allows linking an end to a start
maxAmpDiffSq = (norminv(cutoffProbA1,0,ampDiffStd)).^2;

%calculate the additive constant for each cost as a function of time gap
addConst = log(ampDiffStd) - dispSqR.*log(dispSqTheta) + log(gamma(dispSqR));

%get number of tracks formed by initial linking
numTracks = size(trackedFeatInfo,1);

%get the number of tracks whose starts are to be considered
m = length(indxStart);

%get the number of tracks whose ends are to be considered
n = length(indxEnd);

%get the x,y,z-coordinates and amplitudes of features at the starts of
%their tracks
xCoordStart = zeros(numTracks,1);
yCoordStart = zeros(numTracks,1);
zCoordStart = zeros(numTracks,1);
ampStart = zeros(numTracks,1);
for i=1:numTracks
    xCoordStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*8+1);
    yCoordStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*8+2);
    zCoordStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*8+3);
    ampStart(i) = trackedFeatInfo(i,(trackStartTime(i)-1)*8+4);
end

%get the x,y,z-coordinates and amplitudes of features at the ends of
%their tracks
xCoordEnd = zeros(numTracks,1);
yCoordEnd = zeros(numTracks,1);
zCoordEnd = zeros(numTracks,1);
ampEnd = zeros(numTracks,1);
for i=1:numTracks
    xCoordEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*8+1);
    yCoordEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*8+2);
    zCoordEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*8+3);
    ampEnd(i) = trackedFeatInfo(i,(trackEndTime(i)-1)*8+4);
end

%remove tracks that start at the first time point and get the total number of
%tracks whose starts are to be considered
trackStartTime = trackStartTime(indxStart);
xCoordStart = xCoordStart(indxStart);
yCoordStart = yCoordStart(indxStart);
zCoordStart = zCoordStart(indxStart);
ampStart = ampStart(indxStart);

%remove tracks that end at the last time point and get the total number of
%tracks whose ends are to be considered
trackEndTime = trackEndTime(indxEnd);
xCoordEnd = xCoordEnd(indxEnd);
yCoordEnd = yCoordEnd(indxEnd);
zCoordEnd = zCoordEnd(indxEnd);
ampEnd = ampEnd(indxEnd);

indx1 = []; %row number in cost matrix
indx2 = []; %column number in cost matrix
cost  = []; %cost value

%costs for closing gaps due to features going out-of-focus

for j=1:m %go over all starts
    for i=1:n %go over all ends

        %subtract starting time from ending time
        timeGap = trackStartTime(j) - trackEndTime(i);

        %if the time gap between the end of track i and the start of track
        %j is within the acceptable limits ...
        if timeGap > 1 && timeGap <= timeWindow

            %calculate the square distance between the end
            %point and starting point
            dispSq = (xCoordEnd(i)-xCoordStart(j))^2 + ...
                (yCoordEnd(i)-yCoordStart(j))^2 + ...
                (zCoordEnd(i)-zCoordStart(j))^2;

            %calculate the square difference between the amplitude at
            %the beginning of track j and the end of track i
            ampDiffSq = (ampEnd(i) - ampStart(j))^2;

            %if this is a possible link ...
            if dispSq < maxDispSq(timeGap) && ...
                    ampDiffSq < maxAmpDiffSq(timeGap)

                %assign the cost for this pair of tracks
                indx1 = [indx1; i]; %row number
                indx2 = [indx2; j]; %column number
                cost = [cost; dispSqTheta(timeGap)*dispSq - ...
                    (dispSqR(timeGap)-1)*log(max(dispSq,realmin)) + ...
                    ampDiffSq/2/ampDiffStd(timeGap)^2 + ...
                    addConst(timeGap)]; %cost

            end %(if dispSq < maxDispSq(timeGap) && ampDIffSq < maxAmpDiffSq(timeGap))

        end %(if timeGap > 1 || timeGap <= timeWindow)

    end %(for i=1:n)
end %(for j=1:m)

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting

%if merging and splitting are to be considered ...
if mergeSplit

    %calculate additional additive constants to cost
    halfLog2 = log(2)/2;
    
    %get the maximum squared displacement that allows merging and splitting
    maxDispSqMS = expinv(cutoffProbD2,dispSqR,1./dispSqTheta);

    %get maximum allowed intensity variation when merging or
    %splitting
    maxAmpDiffSqMS = (norminv(cutoffProbA2,0,1.4142*ampDiffStd)).^2;
    
    %costs of merging
    
    %go over all track ending times
    for endTime = min(trackEndTime):max(trackEndTime)

        %find tracks that end at this time point
        tracksToConsider = find(trackEndTime==endTime);
        numTracksToConsider = length(tracksToConsider);

        %first consider the case where the ending track merges with another
        %existing track. This happens when the difference in the intensities
        %of the two merging features is smaller than the change in intensity
        %from one time point to the next

        %get index indicating time of merging
        timeIndx  = endTime*8;

        for j=tracksToConsider' %go over all ends considered
            for i=1:numTracks %go over all tracks

                %get position and amplitude of merged feature
                xCoordMid = trackedFeatInfo(i,timeIndx+1);
                yCoordMid = trackedFeatInfo(i,timeIndx+2);
                zCoordMid = trackedFeatInfo(i,timeIndx+3);
                ampMidT1 = trackedFeatInfo(i,timeIndx+4);

                %get amplitude of feature that merged with the
                %ending track, at time point before merging
                ampMidT0 = trackedFeatInfo(i,timeIndx-3);

                %calculate the square distance between the ending track
                %and the point of merging
                %dispSq = NaN if the track does not exist at the
                %time of merging (prohibiting a link)
                dispSq = (xCoordEnd(j)-xCoordMid)^2 ...
                    + (yCoordEnd(j)-yCoordMid)^2 + (zCoordEnd(j)-zCoordMid)^2;

                %calculate the square difference between the amplitude
                %after merging and the sum of amplitudes before merging
                %ampDiffSq = NaN if the track does not exist at the
                %time of merging or the time point before it (prohibiting a
                %link).
                ampDiffSq = (ampEnd(j) + ampMidT0 - ampMidT1)^2;

                %if this is a possible link ...
                if dispSq < maxDispSqMS(1) && ...
                        ampDiffSq < maxAmpDiffSqMS(1)

                    %increase the "merge index" by one
                    numMerge = numMerge + 1;

                    %save the merging track's number
                    indxMerge = [indxMerge; i];

                    %calculate the cost of merging
                    indx1 = [indx1; j]; %row number
                    indx2 = [indx2; numMerge+m]; %column number

                    cost = [cost; dispSqTheta(1)*dispSq - ...
                        (dispSqR(1)-1)*log(max(dispSq,realmin)) + ...
                        ampDiffSq/4/ampDiffStd(1)^2 + ...
                        addConst(1) + halfLog2]; %cost of merging

                    altCostMerge = [altCostMerge; dispSqR(1) - ...
                        (dispSqR(1)-1)*log(dispSqR(1)/dispSqTheta(1)) + ...
                        (ampMidT0-ampMidT1)^2/4/ampDiffStd(1)^2 + ...
                        addConst(1) + halfLog2]; %cost of not merging

                end %(if dispSq < maxDispSqMS(1) && ampDiffSq < maxAmpDiffSqMS)

            end %(for i=1:numTracks)
        end %(for j=tracksToConsider')

    end %(for endTime = min(trackEndTime):max(trackEndTime))

    %costs of splitting

    %go over all track starting times
    for startTime = min(trackStartTime):max(trackStartTime)
    
        %find tracks that start at this time point
        tracksToConsider = find(trackStartTime==startTime);

        %first consider the case where the starting track splits from another 
        %existing track. This happens when the difference in the intensities
        %of the two merged features is smaller than the change in intensity
        %from one time point to the next

        for j=tracksToConsider' %go over all starts considered
            for i=1:numTracks %go over all tracks

                %get indx indicating time of splitting
                timeIndx  = (trackStartTime(j)-2)*8;

                %get position and amplitude of feature before splitting
                xCoordMid = trackedFeatInfo(i,timeIndx+1);
                yCoordMid = trackedFeatInfo(i,timeIndx+2);
                zCoordMid = trackedFeatInfo(i,timeIndx+3);
                ampMidT1 = trackedFeatInfo(i,timeIndx+4);

                %get amplitude of feature that split at time point
                %after splitting
                ampMidT0 = trackedFeatInfo(i,timeIndx+9);

                %calculate the square distance between the starting track
                %and the point of splitting
                dispSq = (xCoordStart(j)-xCoordMid)^2 ...
                    + (yCoordStart(j)-yCoordMid)^2 + (zCoordStart(j)-zCoordMid)^2;

                %calculate the square difference between the amplitude
                %after merging and the sum of amplitudes before merging
                ampDiffSq = (ampStart(j) + ampMidT0 - ampMidT1)^2;

                %if this is a possible link ...
                if dispSq < maxDispSqMS(1) && ...
                        ampDiffSq < maxAmpDiffSqMS(1)

                    %increase the "split index" by one
                    numSplit = numSplit + 1;

                    %save the splitting track's number
                    indxSplit = [indxSplit; i];

                    %calculate the cost of splitting
                    indx1 = [indx1; numSplit+n]; %row number
                    indx2 = [indx2; j]; %column number

                    cost = [cost; dispSqTheta(1)*dispSq - ...
                        (dispSqR(1)-1)*log(max(dispSq,realmin)) + ...
                        ampDiffSq/4/ampDiffStd(1)^2 + ...
                        addConst(1) + halfLog2]; %cost of splitting

                    altCostSplit = [altCostSplit; dispSqR(i) - ...
                        (dispSqR(1)-1)*log(dispSqR(1)/dispSqTheta(1)) + ...
                        (ampMidT0-ampMidT1)^2/4/ampDiffStd(1)^2 + ...
                        addConst(1) + halfLog2]; %cost of not splitting

                end %(if dispSq < maxDispSqMS(1) && ampDiffSq < maxAmpDiffSqMS)

            end %(for i=1:numTracks)
        end %(for j=tracksToConsider')

    end %(for startTime = min(trackStartTime):max(trackStartTime))
    
end %(if mergeSplit)

%create cost matrix without births and deaths
n1 = n+numSplit;
m1 = m+numMerge;
costMat = sparse(indx1,indx2,cost,n1,m1);

%append cost matrix to allow births and deaths ...

%determine the cost of birth and death
costBD = max(max(max(costMat))+1,1);

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);

%create cost matrix that allows for births and deaths
costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags([costBD*ones(n,1); altCostSplit],0,n1,n1); ... %costs for death
    spdiags([costBD*ones(m,1); altCostMerge],0,m1,m1) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),m1,n1)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);
    
%since the full cost matrix (including the birth and death augmentation) 
%is defined above, the noLnkPrctl variable should not be used for now. 
%determine noLinkCost
if noLnkPrctl ~= -1
    noLinkCost = prctile(cost,noLnkPrctl);
end


%%%%% ~~ the end ~~ %%%%%

