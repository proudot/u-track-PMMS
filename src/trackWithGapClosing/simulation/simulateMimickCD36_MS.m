function [simMPM,tracksSim] = simulateMimickCD36_MS(imSize,numP,lftDist,...
    numF,intVec,mtSpacing,motionParam)
%SIMULATEMIMICKCD36_MS generates tracks that mimick CD36 motion
%
% INPUT 	imSize        : Image size vector [sx,sy]
%           numP          : Average number of points per image.
%           lftDist       : Life time distribution. 
%                           Vector of normalized probability.
%                           The vector is 1-dimensional, as the index
%                           corresponds to the number of frames 
%                           - if e.g. all objects should have the same 
%                           lifetime 10 frames, then lftDist should 
%                           have the form [0 0 0 0 0 0 0 0 0 1].
%           numF          : Number of frames
%           intVec        : Intensity vector [average std]. std refers to the
%                           variation in intensity. In counts (assuming,
%                           for example, a 16-bit camera).
%           mtSpacing     : Mean and standard deviation of spacing between
%                           MTs (in pixels). MTs are distributed along the
%                           x-axis and run parallel to the y-axis.
%           motionParam   : Structure with fields:
%               .diffCoef2D   : Diffusion coefficient of 2D Brownian
%                               motion.
%               .confRad2D    : Confinement radius of 2D Brownian motion.
%               .diffCoef1D   : Diffusion coefficient of 1D Brownian motion
%                               on top of 2D Brownian motion.
%               .fractionLin  : Fraction of tracks that exhibit 1D
%                               diffusion. Skip field or enter [] for
%                               default = 0.5.
%               .probMS       : Row of 4 entries indicating the
%                               probability of having 0, 1, 2 or 3 splits/merges.
%                               The sum of probabilities = 1.
%                               Skip field or enter [] for no merges & splits.
%               .randOrient   : 1 if giving each linear track a random
%                               direction, 0 to give them all the same
%                               direction. Skip field or enter [] for
%                               default = 0;
%
% OUTPUT    simMPM        : Matrix of tracks for Dinah's makeAiryImageFromMPM.
%           tracksFinal   : Tracks in the format of the output of
%                           trackCloseGapsKalman.
%
% Khuloud Jaqaman, September 2007

%%   intialize variables

%get maximum lifetime
lifetimeMax = length(lftDist);

%compute cumulative lifetime distribution function
cumLftDist = zeros(lifetimeMax,1);
cumLftDist(1) = lftDist(1);
for i = 2 : lifetimeMax
    cumLftDist(i) = cumLftDist(i-1) + lftDist(i);
end

% x-length of image
lx = imSize(1);
% y-length of image
ly = imSize(2);

%get motion parameters
diffCoef2D = motionParam.diffCoef2D;
confRad2D = motionParam.confRad2D;
diffCoef1D = motionParam.diffCoef1D;
if isfield(motionParam,'fractionLin') && ~isempty(motionParam.fractionLin)
    fractionLin = motionParam.fractionLin;
else
    fractionLin = 0.5;
end
if isfield(motionParam,'probMS') && ~isempty(motionParam.probMS)
    mergeSplit = 1;
    probMS = motionParam.probMS;
    probMS = [probMS(1) sum(probMS(1:2)) sum(probMS(1:3)) 1];
else
    mergeSplit = 0;
end
if isfield(motionParam,'randOrient') && ~isempty(motionParam.randOrient)
    randOrient = motionParam.randOrient;
else
    randOrient = 0;
end

%%   place MTs in image

%assign average x-coordinate of MTs
if ~isempty(confRad2D)
    mtPosX = (3*confRad2D : mtSpacing(1) : lx-3*confRad2D)';
else
    mtPosX = (1 : mtSpacing(1) : lx-1)';
end

%get number of MTs
numMT = length(mtPosX);

%perturb MT x-coordinates so that distribution is not uniform
mtPosX = mtPosX + randn(numMT,1)*mtSpacing(2);


%%   determine number of iterations based on number of objects and frames

% expectancy value for lifetime
ex_lft = sum((1:lifetimeMax).*shiftdim(lftDist)');

% the necessary number of simulated tracks to reach the required specified
% density of objects per image is approximately
numTracks = round(numP*numF/ex_lft); 


%%   create objects with specified lifetime, initial position & intensity, motion type and merge/split events

%initialize tracksSim
tracksSim = repmat(struct('tracksCoordAmpCG',[],'tracksFeatIndxCG',[],'seqOfEvents',[]),numTracks,1);

%assign track start times, end times and lifetimes
numTracksTmp = 0;
startframe = [];
lifetime = [];
endframe = [];
while numTracksTmp < numTracks

    %randomly choose 10*numTracks track starting frames
    %allow the search to go back "maximum lifetime" frames before start of the movie
    startframeTmp = round((lifetimeMax+numF-2)*rand(10*numTracks,1)) - lifetimeMax + 2;

    %assign track lifetimes based on the input lifetime distribution
    randVar = rand(10*numTracks,1);
    randVar = randVar(randVar<max(cumLftDist));
    numAttempts = length(randVar);
    lifetimeTmp = zeros(numAttempts,1);
    for iTrack = 1 : numAttempts
        lifetimeTmp(iTrack) = find(cumLftDist>=randVar(iTrack),1,'first');
    end

    %hence calculate ending frame
    startframeTmp = startframeTmp(1:numAttempts);
    endframeTmp = startframeTmp + lifetimeTmp - 1;

    %retain only tracks that end after frame 1
    indxKeep = find(endframeTmp >= 1);
    numKeep = length(indxKeep);
    startframeTmp = startframeTmp(indxKeep);
    lifetimeTmp = lifetimeTmp(indxKeep);
    endframeTmp = endframeTmp(indxKeep);

    %retain only the first numTracks tracks if there are that many
    startframe = [startframe; startframeTmp(1:min(numKeep,numTracks))];
    lifetime = [lifetime; lifetimeTmp(1:min(numKeep,numTracks))];
    endframe = [endframe; endframeTmp(1:min(numKeep,numTracks))];
    
    %get number of tracks generated
    numTracksTmp = length(startframe);
    
end %(while numTracksTmp < numTracks)

%truncate start and end times (and consequently lifetimes) to be within the
%movie
vis_startframe = max([startframe ones(numTracks,1)],[],2);
vis_endframe   = min([endframe numF*ones(numTracks,1)],[],2);
vis_lifetime   = vis_endframe - vis_startframe + 1;

%assign initial positions
posX = mtPosX(ceil(rand(numTracks,1)*numMT)); %x-coordinate - on one of the MTs
if ~isempty(confRad2D)
    posY = 3*confRad2D + rand(numTracks,1) * (ly - 6*confRad2D); %y-coordinate - anywhere on an MT
else
    posY = 1 + rand(numTracks,1) * (ly - 2);
end

%assign initial intensities
startInt = intVec(1) + intVec(2)*randn(numTracks,1);

%assign motion types (0: Brownian, 1: Brownian + linear)
mType = rand(numTracks,1) <= fractionLin;

%go over tracks and assign them a sequence of events, i.e. start times,
%end times and merge/split times
%follow the convention used in trackCloseGapsKalman in documenting the
%sequence of events
%for now, every split is immediately followed by a merge
randVar = rand(numTracks,1);
for iTrack = 1 : numTracks

    %first store the start and end time of the main track
    seqOfEvents = [vis_startframe(iTrack) 1 1 NaN; vis_endframe(iTrack) 2 1 NaN];

    %then include merges and splits
    if mergeSplit && vis_lifetime(iTrack) > 10

        %based on the chosen random numbers, decide the number of merges
        %and splits
        numMS = (find(probMS>=randVar(iTrack),1,'first')-1);

        %if there are merges and splits
        if numMS > 0

            %assign times when neighbor appears
            timeMS = randsample((vis_startframe(iTrack)+2:2:vis_endframe(iTrack)-2)',numMS);

            %add merges and splits to the sequence of events
            seqOfEvents = [seqOfEvents; ...
                [timeMS(1:numMS) ones(numMS,1) (1:numMS)'+1 ones(numMS,1)]; ... %all the splits
                [timeMS(1:numMS)+1 2*ones(numMS,1) (1:numMS)'+1 ones(numMS,1)]]; %all the merges

            %sort sequence of events in ascending order of time
            [dummy,indx] = sort(seqOfEvents(:,1));
            seqOfEvents = seqOfEvents(indx,:);

        end
        
    end %(if mergeSplit && vis_lifetime(iTrack) > 10)

    %store sequence of events in tracksSim
    tracksSim(iTrack).seqOfEvents = seqOfEvents;

end %(for iTrack = 1 : numTracks)

%%   create trajectories for all objects in the list, and store them in tracksSim

%go over all objects ...
for iTrack = 1 : numTracks
    
    %get sequence of events for this object
    seqOfEvents = tracksSim(iTrack).seqOfEvents;
    
    %current number of needed frames 
    cnf = vis_lifetime(iTrack);
    
    %start position of this object
    xystart = [posX(iTrack) posY(iTrack)];
    
    %generate main track
    if cnf > 1
        if ~isempty(confRad2D)
            xyvecTraj = brownianMotion(2,diffCoef2D,cnf-1,0.1,1,confRad2D); %Brownian part
        else
            xyvecTraj = brownianMotion(2,diffCoef2D,cnf-1,0.1);
        end
        if mType(iTrack) == 1 %linear part
            trackLin = brownianMotion(1,diffCoef1D,cnf-1,0.1);
            %             trackLinDiff = (1+0.1*(2*rand((cnf-1)*10,1)-1)).*sqrt(2*diffCoef1D*0.1)...
            %                 .*sign(randn((cnf-1)*10,1));
            %             trackLin = zeros((cnf-1)*10+1,1);
            %             for i=1:length(trackLin)-1
            %                 trackLin(i+1) = trackLin(i) + trackLinDiff(i);
            %             end
            if randOrient
                orientAngle = rand(1)*2*pi;
                trackLinXY = [trackLin*cos(orientAngle) trackLin*sin(orientAngle)];
                xyvecTraj = xyvecTraj + trackLinXY;
            else
                xyvecTraj(:,2) = xyvecTraj(:,2) + trackLin;
            end
        end
        xyvecTraj = xyvecTraj(1:10:end,:);
    else
        xyvecTraj = zeros(1,2);
    end
    xyvecTraj = xyvecTraj + repmat(xystart,cnf,1); %add initial position
    
    %assign track intensity
    intVecTraj = [startInt(iTrack); intVec(1)+intVec(2)*randn(cnf-1,1)];
    
    %don't allow nonpositive intensities
    intVecTraj(intVecTraj<=0) = 1;

    %if there are possible merges and splits
    if mergeSplit
        
        %get splitting times and number of merges/splits
        indxMS = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
        numMS = length(indxMS);
        timeMS = seqOfEvents(indxMS,1);
        
        %assign positions
        posMS = NaN(numMS,2);
        for iMS = 1 : numMS
            mainTrackPos = xyvecTraj(timeMS(iMS)-vis_startframe(iTrack):...
                timeMS(iMS)-vis_startframe(iTrack)+1,:);
            posMS(iMS,:) = mainTrackPos(2,:) + (2*randn(1,2)-1) .* ...
                (mainTrackPos(2,:)-mainTrackPos(1,:))/50;
        end

        %add up intensities to account for merging and splitting
        intVecMS = intVecTraj(timeMS-vis_startframe(iTrack)+1) .* (rand(numMS,...
            1)*0.2+0.4);
        intVecTraj(timeMS-vis_startframe(iTrack)+1) = intVecTraj(...
            timeMS-vis_startframe(iTrack)+1) .* (rand(numMS,1)*0.2+0.4);

        %allocate memory for tracksCoordAmpCG
        tracksCoordAmpCG = NaN(numMS+1,8*cnf);
        
        %info of main track
        tracksCoordAmpCG(1,1:8:end) = xyvecTraj(:,1);
        tracksCoordAmpCG(1,2:8:end) = xyvecTraj(:,2);
        tracksCoordAmpCG(1,4:8:end) = intVecTraj;
        
        %info for merging/splitting branches
        for iMS = 1 : numMS
            iSegment = seqOfEvents(indxMS(iMS),3);
            tracksCoordAmpCG(iSegment,(timeMS(iMS)-vis_startframe(iTrack))*8+1:...
                (timeMS(iMS)-vis_startframe(iTrack))*8+2) = posMS(iMS,:);
            tracksCoordAmpCG(iSegment,(timeMS(iMS)-vis_startframe(iTrack))*8+4:...
                (timeMS(iMS)-vis_startframe(iTrack))*8+4) = intVecMS(iMS);
        end

    else %no merges and splits

        tracksCoordAmpCG = NaN(1,8*cnf);
        tracksCoordAmpCG(1,1:8:end) = xyvecTraj(:,1);
        tracksCoordAmpCG(1,2:8:end) = xyvecTraj(:,2);
        tracksCoordAmpCG(1,4:8:end) = intVecTraj;
        
    end

    %savetrack in structure
    tracksSim(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;

end

%%   make simMPM out of tracksFinal

%allocate memory for simMPM
simMPM = zeros(numTracks,8*numF);

%get number of segments making each track
numSegments = zeros(numTracks,1);
for i = 1 : numTracks
    numSegments(i) = size(tracksSim(i).tracksCoordAmpCG,1);
end

%locate the row of the first track of each compound track in the
%big matrix of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%put all tracks together in a matrix
for i = 1 : numTracks
    startTime = tracksSim(i).seqOfEvents(1,1);
    endTime   = tracksSim(i).seqOfEvents(end,1);
    simMPM(trackStartRow(i):trackStartRow(i)+...
        numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
        tracksSim(i).tracksCoordAmpCG;
end

%remove extra columns
dummy = simMPM;
clear simMPM
simMPM = zeros(size(dummy,1),3*numF);
simMPM(:,1:3:end) = dummy(:,1:8:end);
simMPM(:,2:3:end) = dummy(:,2:8:end);
simMPM(:,3:3:end) = dummy(:,4:8:end);
simMPM(isnan(simMPM)) = 0;

%% shift positions to remove any negative coordinates

%get minimum x and y coordinates
minX = min(min(simMPM(:,1:3:end)));
minY = min(min(simMPM(:,2:3:end)));

%shift if minimum is nonpositive
if minX <= 0
    shiftX = abs(minX) + 1;
    simMPM(:,1:3:end) = simMPM(:,1:3:end) + shiftX;
    for iTrack = 1 : numTracks
        tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end) = ...
            tracksSim(iTrack).tracksCoordAmpCG(:,1:8:end) + shiftX;
    end
end

if minY <= 0
    shiftY = abs(minY) + 1;
    simMPM(:,2:3:end) = simMPM(:,2:3:end) + shiftY;
    for iTrack = 1 : numTracks
        tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end) = ...
            tracksSim(iTrack).tracksCoordAmpCG(:,2:8:end) + shiftY;
    end
end

