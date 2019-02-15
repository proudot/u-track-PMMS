function [simMPM,tracksSim,trackDiffCoef] = simulateRandomPlusLinearMotionHetero(imSize,numP,...
    lftDist,numF,intVec,motionParam)
%SIMULATERANDOMPLUSLINEARMOTION generates tracks that can exhibit unconfined, confined and directed motion
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
%           motionParam   : Structure with fields:
%               .diffCoef2D   : 2-element row vector indicating range of
%                               diffusion coefficients for 2D Brownian
%                               motion.
%               .confRad2D    : 2-element row vector indicating range of
%                               confinement radii for 2D confined motion.
%               .speed1D      : 2-element row vector inficating range of
%                               speeds for directed motion. Note that
%                               direction can be anything and will be
%                               chosen randomly.
%               .probVelSwitch: Probability of switching velocities for
%                               directed particles (magnitude and
%                               direction).
%               .fracType     : Fraction of tracks exhibiting unconfined,
%                               confined and directed motion.
%               .probMS       : Row of 4 entries indicating the
%                               probability of having 0, 1, 2 or 3 splits/merges.
%                               The sum of probabilities = 1.
%                               Skip field or enter [] for no merges & splits.
%
% OUTPUT    simMPM        : Matrix of tracks for Dinah's makeAiryImageFromMPM.
%           tracksFinal   : Tracks in the format of the output of
%                           trackCloseGapsKalman.
%
% Khuloud Jaqaman, September 2011

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
speed1D = motionParam.speed1D;
probVelSwitch = motionParam.probVelSwitch;
probTypeSwitch = motionParam.probTypeSwitch;
fracType = motionParam.fracType;
if isfield(motionParam,'probMS') && ~isempty(motionParam.probMS)
    mergeSplit = 1;
    probMS = motionParam.probMS;
    probMS = [probMS(1) sum(probMS(1:2)) sum(probMS(1:3)) 1];
else
    mergeSplit = 0;
end

%%   determine number of iterations based on number of objects and frames

% expectancy value for lifetime
ex_lft = sum((1:lifetimeMax).*shiftdim(lftDist)');

% the necessary number of simulated tracks to reach the required specified
% density of objects per image is approximately
numTracks = round(numP*numF/ex_lft); 


%%   create objects with specified lifetime, initial position & intensity, motion parameters and merge/split events

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
    
    %make sure that starts and ends are within the movie time
    startframeTmp = max(startframeTmp,1);
    endframeTmp = min(endframeTmp,numF);
    lifetimeTmp = endframeTmp - startframeTmp + 1;

    %retain only tracks that end after frame 1
    %also retain only tracks that are at least 2 frames long
    indxKeep = find(endframeTmp >= 1 & lifetimeTmp >= 2);
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
posX = 1 + rand(numTracks,1) * (lx - 2);
posY = 1 + rand(numTracks,1) * (ly - 2);

%assign initial intensities
startInt = intVec(1) + intVec(2)*randn(numTracks,1);

%assign motion types (0: Brownian, 1: Brownian + linear)
randNumber = rand(numTracks,1);
mType = zeros(numTracks,1); %start with unconfined
cumFrac=cumsum(fracType);
for i = (size(cumFrac,2)):-1:1
    mType(randNumber<cumFrac(i)) = i;
end

%assign motion parameters
trackDiffCoef = rand(numTracks,1)*diff(diffCoef2D) + diffCoef2D(1);
trackConfRad = rand(numTracks,1)*diff(confRad2D) + confRad2D(1);
timeStep=0.1;



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

            %get time points available for splitting
            availableTimes = (vis_startframe(iTrack)+6:vis_endframe(iTrack)-14)';
            
            if length(availableTimes) >= numMS
                
                %assign splitting times
                if length(availableTimes)==1
                    timeSplit = availableTimes;
                else
                    timeSplit = randsample(availableTimes,numMS);
                end
                
                %assign merging times
                timeMerge = timeSplit + 8;
%                 timeMerge = timeSplit + round(rand(numMS,1)*3+5);
                
                %if any merges take place after the end of the track, remove
                %both split and merge
                indxGood = find(timeMerge < vis_endframe(iTrack));
                timeSplit = timeSplit(indxGood);
                timeMerge = timeMerge(indxGood);
                numMS = length(indxGood);
                
                if numMS > 0
                    
                    %add merges and splits to the sequence of events
                    seqOfEvents = [seqOfEvents; ...
                        [timeSplit ones(numMS,1) (1:numMS)'+1 ones(numMS,1)]; ... %all the splits
                        [timeMerge 2*ones(numMS,1) (1:numMS)'+1 ones(numMS,1)]]; %all the merges
                    
                    %sort sequence of events in ascending order of time
                    [dummy,indx] = sort(seqOfEvents(:,1));
                    seqOfEvents = seqOfEvents(indx,:);
                    
                end
                
            end

        end
        
    end %(if mergeSplit && vis_lifetime(iTrack) > 10)

    %store sequence of events in tracksSim
    tracksSim(iTrack).seqOfEvents = seqOfEvents;

end %(for iTrack = 1 : numTracks)

%%   create trajectories for all objects in the list, and store them in tracksSim

%initialize vector to store indices of tracks to remove in the end
trackGood = ones(numTracks,1);

%go over all objects ...
for iTrack = 1 : numTracks
    
    %get sequence of events for this object
    seqOfEvents = tracksSim(iTrack).seqOfEvents;
    
    %current number of needed frames 
    cnf = vis_lifetime(iTrack);
    
    %start position of this object
    xystart = [posX(iTrack) posY(iTrack)];
    
    %track type
    trackType = mType(iTrack);
    
    %generate main track
    if cnf > 1
        
        switch trackType
            
          case 1 %unconfined
            
            %get trajectory
            xyvecTraj = brownianMotion(2,trackDiffCoef(iTrack),cnf-1,timeStep);
            typeVecTraj=1*ones(1,size(xyvecTraj,1));
            
          case 2 %confined
            
            %get trajectory
            xyvecTraj = brownianMotion(2,trackDiffCoef(iTrack),cnf-1,timeStep,1,trackConfRad(iTrack));
            typeVecTraj=2*ones(1,size(xyvecTraj,1));
            
          case 3 %directed
            
            %determine time points at which a new velocity is assigned
            timeVelSwitch = [1; 1+find(rand(cnf-2,1)<probVelSwitch)];
            
            %determine number of segments and length of each segment
            numSeg = length(timeVelSwitch);
            lengthSeg = diff([timeVelSwitch; cnf]) + 1;
            
            %assign velocities
            trackVel = [rand(numSeg,1)*diff(speed1D)+speed1D(1) rand(numSeg,1)*360];
            
            %get trajectory, in segments
            xyvecTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(1),timeStep,[],[],trackVel(1,:));
            typeVecTraj=3*ones(1,size(xyvecTraj,1));
            for iSeg = 2 : numSeg
                tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep,[],[],trackVel(iSeg,:));
                tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
                tmpTypeTraj=3*ones(1,size(tmpTraj,1));
                xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
                typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
            end
            
          case 4 %Heterogenous evrything goes
            
            %determine time points at which a new velocity is assigned
            timeVelSwitch = [1; 1+find(rand(cnf-2,1)<probVelSwitch)];
            
            %determine number of segments and length of each segment
            numSeg = length(timeVelSwitch);
            lengthSeg = diff([timeVelSwitch; cnf]) + 1;
            
            %assign velocities
            trackVel = [rand(numSeg,1)*diff(speed1D)+speed1D(1) rand(numSeg,1)*360];
            
            %start with jiggling motion
            xyvecTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(1),timeStep,1,trackConfRad(iTrack));
            typeVecTraj=2*ones(1,size(xyvecTraj,1));
            motion_type=1;
            for iSeg = 2 : numSeg
                switch motion_type
                  case 1 %unconfined
                    tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep);
                    tmpTypeTraj=1*ones(1,size(tmpTraj,1));
                  case 2 %confined
                    tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep,1,trackConfRad(iTrack));
                    tmpTypeTraj=2*ones(1,size(tmpTraj,1));
                  case 3 %directed 
                    tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep,[],[],trackVel(iSeg,:));
                    tmpTypeTraj=3*ones(1,size(tmpTraj,1));
                end 
                motion_type=mod(motion_type,3)+1;
                tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
                xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
                typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
            end     
            
          case 5 %Heterogenous unconfined+directed
            
            %determine time points at which a new velocity is assigned
            timeVelSwitch = [1; 1+find(rand(cnf-2,1)<probTypeSwitch)];
            
            %determine number of segments and length of each segment
            numSeg = length(timeVelSwitch);
            lengthSeg = diff([timeVelSwitch; cnf]) + 1;
            
            %assign velocities

            
            %start with jiggling motion (confined)
            xyvecTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(1),timeStep,1,trackConfRad(iTrack));
            typeVecTraj=2*ones(1,size(xyvecTraj,1));
            motion_type=1;
            for iSeg = 2 : numSeg
                switch motion_type
                  case 1 %unconfined
                    tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep);
                    tmpTypeTraj=1*ones(1,size(tmpTraj,1));
                  case 2 %directed 
                    %determine time points at which a new velocity is assigned
                    timeVelSwitch = [1; 1+find(rand(lengthSeg(iSeg)-2,1)<probVelSwitch)];
                    %determine number of segments and length of each segment
                    numSubSeg = length(timeVelSwitch);
                    lengthSubSeg = diff([timeVelSwitch;lengthSeg(iSeg)]) + 1;
                    trackVel = [rand(numSubSeg,1)*diff(speed1D)+speed1D(1) rand(numSubSeg,1)*360];
                    
                    xyvecSubTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(1),timeStep,[],[],trackVel(1,:));
                    for iSubSeg = 2:numSubSeg
                        tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(iSubSeg),timeStep,[],[],trackVel(iSubSeg,:));
                        tmpTraj = tmpTraj + repmat(xyvecSubTraj(end,:),size(tmpTraj,1),1);
                        xyvecSubTraj = [xyvecSubTraj; tmpTraj(2:end,:)];
                    end 
                    tmpTraj=xyvecSubTraj;
                    tmpTypeTraj=3*ones(1,size(tmpTraj,1));
                end 
                motion_type=mod(motion_type,2)+1;
                tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
                xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
                typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
            end    
          case 6 %Heterogenous confined+directed
            
            %determine time points at which a new velocity is assigned
            timeVelSwitch = [1; 1+find(rand(cnf-2,1)<probTypeSwitch)];
            
            %determine number of segments and length of each segment
            numSeg = length(timeVelSwitch);
            lengthSeg = diff([timeVelSwitch; cnf]) + 1;
            
            %assign velocities
            
            %start with jiggling motion
            xyvecTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(1),timeStep,1,trackConfRad(iTrack));
            typeVecTraj=2*ones(1,size(xyvecTraj,1));
            motion_type=2;
            for iSeg = 2 : numSeg
                switch motion_type
                  case 1 %confined
                    tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep,1,trackConfRad(iTrack));
                    tmpTypeTraj=2*ones(1,size(tmpTraj,1));
                  case 2 %directed 
                         %determine time points at which a new velocity is assigned
                    timeVelSwitch = [1; 1+find(rand(lengthSeg(iSeg)-2,1)<probVelSwitch)];
                    %determine number of segments and length of each segment
                    numSubSeg = length(timeVelSwitch);
                    lengthSubSeg = diff([timeVelSwitch;lengthSeg(iSeg)]) + 1;
                    trackVel = [rand(numSubSeg,1)*diff(speed1D)+speed1D(1) rand(numSubSeg,1)*360];
                    
                    xyvecSubTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(1),timeStep,[],[],trackVel(1,:));
                    for iSubSeg = 2:numSubSeg
                        tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(iSubSeg),timeStep,[],[],trackVel(iSubSeg,:));
                        tmpTraj = tmpTraj + repmat(xyvecSubTraj(end,:),size(tmpTraj,1),1);
                        xyvecSubTraj = [xyvecSubTraj; tmpTraj(2:end,:)];
                    end 
                    tmpTraj=xyvecSubTraj;
                    tmpTypeTraj=3*ones(1,size(tmpTraj,1));
                end 
                motion_type=mod(motion_type,2)+1;
                tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
                xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
                typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
            end    
          case 7 %Heterogenous directed +  confined motion of half size.
                 %determine time points at which a new velocity is assigned
            timeVelSwitch = [1; 1+find(rand(cnf-2,1)<probTypeSwitch)];
            
            %determine number of segments and length of each segment
            numSeg = length(timeVelSwitch);
            lengthSeg = diff([timeVelSwitch; cnf]) + 1;
            
            %start with jiggling motion
            xyvecTraj = [];
            typeVecTraj= [];
            motion_type=2;
            for iSeg = 1 : numSeg
                switch motion_type
                  case 1 %confined
                    tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSeg(iSeg),timeStep,1,trackConfRad(iTrack));
                    tmpTypeTraj=2*ones(1,size(tmpTraj,1));
                  case 2 %directed 
                         %determine time points at which a new velocity is assigned
                    timeVelSwitch = [1; 1+find(rand(lengthSeg(iSeg)-2,1)<probVelSwitch)];
                    %determine number of segments and length of each segment
                    numSubSeg = length(timeVelSwitch);
                    lengthSubSeg = diff([timeVelSwitch;lengthSeg(iSeg)]) + 1;
                    trackVel = [rand(numSubSeg,1)*diff(speed1D)+speed1D(1) rand(numSubSeg,1)*360];
                    
                    xyvecSubTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(1),timeStep,[],[],trackVel(1,:));
                    for iSubSeg = 2:numSubSeg
                        tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(iSubSeg),timeStep,[],[],trackVel(iSubSeg,:));
                        tmpTraj = tmpTraj + repmat(xyvecSubTraj(end,:),size(tmpTraj,1),1);
                        xyvecSubTraj = [xyvecSubTraj; tmpTraj(2:end,:)];
                    end 
                    tmpTraj=xyvecSubTraj;
                    tmpTypeTraj=3*ones(1,size(tmpTraj,1));
                end 
                motion_type=mod(motion_type,2)+1;
                if(~isempty(xyvecTraj))
                    tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
                    xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
                    typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
                else
                    xyvecTraj=tmpTraj;
                    typeVecTraj=tmpTypeTraj;
                end                 
            end               
% $$$             %determine time points at which a new velocity is assigned
% $$$             timeVelSwitch = [1; 1+find(rand(cnf-2,1)<probTypeSwitch)];
% $$$             
% $$$             %determine number of segments and length of each segment
% $$$             numSeg = length(timeVelSwitch);
% $$$             lengthSeg = diff([timeVelSwitch; cnf]) + 1;
% $$$ 
% $$$             %assign velocities
% $$$             xyvecTraj=[];
% $$$             typeVecTraj=[];
% $$$             motion_type=2;
% $$$             for iSeg = 1 : numSeg
% $$$                 timeVelSwitch = [1; 1+find(rand(floor(lengthSeg(iSeg)/2)-2,1)<probVelSwitch)];
% $$$                 %determine number of segments and length of each segment
% $$$                 numSubSeg = length(timeVelSwitch);
% $$$                 lengthSubSeg = diff([timeVelSwitch;floor(lengthSeg(iSeg)/2)]) + 1;
% $$$                 trackVel = [rand(numSubSeg,1)*diff(speed1D)+speed1D(1) rand(numSubSeg,1)*360];
% $$$                 
% $$$                 xyvecSubTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(1),timeStep,[],[],trackVel(1,:));
% $$$                 for iSubSeg = 2:numSubSeg
% $$$                     tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),lengthSubSeg(iSubSeg),timeStep,[],[],trackVel(iSubSeg,:));
% $$$                     tmpTraj = tmpTraj + repmat(xyvecSubTraj(end,:),size(tmpTraj,1),1);
% $$$                     xyvecSubTraj = [xyvecSubTraj; tmpTraj(2:end,:)];
% $$$                 end 
% $$$                 tmpTraj=xyvecSubTraj;
% $$$                 tmpTypeTraj=3*ones(1,size(tmpTraj,1));
% $$$ 
% $$$                 if(~isempty(xyvecTraj))
% $$$                     tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
% $$$                     xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
% $$$                     typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
% $$$                 else
% $$$                     xyvecTraj=tmpTraj;
% $$$                     typeVecTraj=tmpTypeTraj;
% $$$                 end                 
% $$$ 
% $$$                 tmpTraj = brownianMotion(2,trackDiffCoef(iTrack),round(lengthSeg(iSeg)/2),timeStep,1,trackConfRad(iTrack));
% $$$                 tmpTypeTraj=2*ones(1,size(tmpTraj,1));
% $$$ 
% $$$                 tmpTraj = tmpTraj + repmat(xyvecTraj(end,:),size(tmpTraj,1),1);
% $$$                 xyvecTraj = [xyvecTraj; tmpTraj(2:end,:)];
% $$$                 typeVecTraj = [typeVecTraj(1:(end-1)) tmpTypeTraj];
% $$$             end   
        end %(switch trackType)
        
        %sample trajectory at time points equivalent to movie frames
        xyvecTraj = xyvecTraj(1:10:end,:);
        xyvecTraj = xyvecTraj(1:cnf,:);

        typeVecTraj= typeVecTraj(1:10:end);
        typeVecTraj = typeVecTraj(1:cnf);

    else
        
        xyvecTraj = zeros(1,2);
        
    end %(if cnf > 1)
    
    %add initial position
    xyvecTraj = xyvecTraj + repmat(xystart,cnf,1);
    
    %find first point outside of the image area
    indxBad = find( xyvecTraj(:,1)>lx | xyvecTraj(:,1)<1 | ...
        xyvecTraj(:,2)>ly | xyvecTraj(:,2)<1,1,'first');
    
    %if there are points outside the image area
    if ~isempty(indxBad)
        
        %remove them
        xyvecTraj = xyvecTraj(1:indxBad-1,:);
        typeVecTraj = typeVecTraj(1:indxBad-1);

        %update lifetime
        cnf = size(xyvecTraj,1);
        
        %if lifetime is still larger than one
        if cnf > 1
            
            %update lifetime and end frame
            vis_lifetime(iTrack) = cnf;
            vis_endframe(iTrack) = vis_startframe(iTrack) + cnf - 1;
            
            %update sequence of events to remove any events happening after
            %the end of the track
            indxBad = find(seqOfEvents(:,1)>vis_endframe(iTrack));
            seqOfEvents(indxBad,1) = vis_endframe(iTrack);
            seqOfEvents(indxBad,4) = NaN;
            indxBadStart = seqOfEvents(seqOfEvents(:,1)==vis_endframe(iTrack)&...
                seqOfEvents(:,2)==1,3);
            for i = indxBadStart'
                seqOfEvents(seqOfEvents(:,3)==i,:) = [];
            end
            tracksSim(iTrack).seqOfEvents = seqOfEvents;
            
        else
            
            trackGood(iTrack) = 0;
            seqOfEvents = seqOfEvents(seqOfEvents(:,3)==1,:);
            
        end %(if cnf > 1)
        
    end %(if ~isempty(indxBad))
    
    %assign track intensity
    intVecTraj = [startInt(iTrack); intVec(1)+intVec(2)*randn(cnf-1,1)];
    
    %don't allow nonpositive intensities
    intVecTraj(intVecTraj<=0) = 1;

    %if there are possible merges and splits
    if mergeSplit
        
        %get splitting times and number of splits
        indxSplits = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
        numSplits = length(indxSplits);
        timeSplits = seqOfEvents(indxSplits,1);
        
        %allocate memory for tracksCoordAmpCG
        tracksCoordAmpCG = NaN(numSplits+1,8*cnf);
        
        %info of main track
        %intensity info will be modified if there are splits and merges
        tracksCoordAmpCG(1,1:8:end) = xyvecTraj(:,1);
        tracksCoordAmpCG(1,2:8:end) = xyvecTraj(:,2);
        tracksCoordAmpCG(1,4:8:end) = intVecTraj;
        
        %info of splits/merges
        for iSplit = 1 : numSplits
            
            iSegment = seqOfEvents(indxSplits(iSplit),3);
            
            %find merging/ending time
            timeEnd = seqOfEvents(seqOfEvents(:,3)==iSegment & seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)),1) - 1;
            if isempty(timeEnd)
                timeEnd = vis_endframe(iTrack);
            end
            
            %store information
            tracksCoordAmpCG(1,(timeSplits(iSplit)-vis_startframe(iTrack))*8+4:8:...
                (timeEnd-vis_startframe(iTrack)+1)*8) = ...
                tracksCoordAmpCG(1,(timeSplits(iSplit)-vis_startframe(iTrack))*8+4:8:...
                (timeEnd-vis_startframe(iTrack)+1)*8) / 2;
            tracksCoordAmpCG(iSegment,(timeSplits(iSplit)-vis_startframe(iTrack))*8+1:...
                (timeEnd-vis_startframe(iTrack)+1)*8) = ...
                tracksCoordAmpCG(1,(timeSplits(iSplit)-vis_startframe(iTrack))*8+1:...
                (timeEnd-vis_startframe(iTrack)+1)*8);
        end
        
    else %no merges and splits
        
        tracksCoordAmpCG = NaN(1,8*cnf);
        tracksCoordAmpCG(1,1:8:end) = xyvecTraj(:,1);
        tracksCoordAmpCG(1,2:8:end) = xyvecTraj(:,2);
        tracksCoordAmpCG(1,4:8:end) = intVecTraj;
        
    end

    %save track in structure
    tracksSim(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
    tracksSim(iTrack).tracksSegType=typeVecTraj;
end

%keep only good tracks
tracksSim = tracksSim(trackGood==1);
numTracks = length(tracksSim);

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

