function simMPM = simulateMimickCD36(imSize,numP,lftDist,numF,intVec,...
    mtSpacing,motionParam)
%SIMULATEMIMICKCD36 generates tracks that mimick CD36 motion
%
% INPUT 	imSize        : Image size vector [sx,sy]
%           numP          : Average number of points per image.
%           lftDist       : Lifetime distribution vector (normalized 
%                           probability). the vector is 1-dimensional, as the
%                           length position automatically corresponds to the 
%                           number of frames - if e.g. all objects should
%                           have the same lifetime 10 frames, then lftDist
%                           should have the form [0 0 0 0 0 0 0 0 0 1]           
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
%                               diffusion.
%
% OUTPUT   MPM file (containing x,y, and intensity), which can later be
%          transformed into a simulated movies with the function makeAiryImageFromMPM
%
% Khuloud & Dinah, September 2007

%%   intialize variables

% maximum lifetime
pl = length(lftDist);

% x-length of image
lx = imSize(1);
% y-length of image
ly = imSize(2);

%get motion parameters
diffCoef2D = motionParam.diffCoef2D;
confRad2D = motionParam.confRad2D;
diffCoef1D = motionParam.diffCoef1D;
fractionLin = motionParam.fractionLin;

%%   place MTs in image

%assign average x-coordinate of MTs
mtPosX = (1.5*confRad2D : mtSpacing(1) : lx-1.5*confRad2D)';

%get number of MTs
numMT = length(mtPosX);

%perturb MT x-coordinates so that distribution is not uniform
mtPosX = mtPosX + randn(numMT,1)*mtSpacing(2);


%%   determine number of iterations based on number of objects and frames

% expectancy value for lifetime
ex_lft = sum((1:pl).*shiftdim(lftDist)');

% the necessary number of simulated objects to reach the required specified
% density of objects per image is approximately
ni = round(numP*numF/ex_lft); 


%%   create objects with specified lifetime, initial position & intensity and motion type

%pre-allocate memory
objectList = repmat(struct('startframe',[],'lft',[],'startpos',[],'startint',[],'mType',[]),ni,1);
allLFT = zeros(ni,1);

% loop over ni
i=1; %particle index
while i <= ni
    
    % randomly select object's starting frame, allow the search to go back
    % 100 frames before start of the movie
    startframe = round((100+numF)*rand(1))-100;
    
    % randomly select lifetime of this object, based on the specified
    % lifetime distribution
    ptest=1; 
    pdist=0;
    while pdist<ptest
        % choose random lifetime between 1 and pl
        plft  = randsample(pl,1);
        % probability for this lifetime
        pdist = lftDist(plft);
        % random number
        ptest = rand(1);
        % passtest is pdist<ptest
    end
    
    % the lifetime randomly assigned to this object is plft; if the object
    % with this lifetime is still visible in frame 1, then enter this 
    % object into the object list and increase the counter
    endframe = (startframe+plft-1);
    
    if endframe>0
        
        % note: in the object list, we only consider the visible lifetime 
        % and the visible startpoint in the movie (not e.g. the real
        % startpoint which may be before the first frame), so we have to 
        % cut away all the frames before the start or after the end of the 
        % movie
        vis_startframe = max(startframe,1);
        vis_endframe   = min(numF,endframe);
        vis_lifetime   = vis_endframe - vis_startframe + 1;
        
        % enter values into object list
        objectList(i).startframe = vis_startframe;
        objectList(i).lft = vis_lifetime;
        
        % assign this object a random initial position
        px = mtPosX(ceil(rand(1)*numMT)); %x-coordinate - on one of the MTs
        py = rand(1)*ly; %y-coordinate - anywhere on an MT
        objectList(i).startpos = [px py];
        
        % assign this object a random start intensity
        objectList(i).startint = intVec(1)+intVec(2)*randn(1);
        
        %assign this object a motion type (0: Brownian, 1: Brownian + linear)
        objectList(i).mType = rand(1) <= fractionLin;

        % enter the actual lifetime into the allLFT vector, which allows us
        % to check the simulated lifetime distribution later on if
        % necessary
        allLFT(i) = plft;
        
        %add one to particle index
        i = i + 1;
        
    end %(if endframe>0)
    
end %(while i<=ni)

% % % CHECKPOINT: if desired, compare the resulting lifetime distribution to
% % % the original probability density at this point
% % allLFThist = hist(allLFT,(1:pl));
% % allLFTnorm = allLFThist/sum(allLFThist);
% % % comment/uncomment the following paragraph to display the simulated
% % % lifetime distribution
% % figure;
% % bar(allLFTnorm);
% % hold on;
% % plot(lftDist,'r.-');


%%   create trajectories for all objects in the list, and enter them into the MPM files

% initialize mpm file 
simMPM = zeros(ni,3*numF);

%go over all objects ...
for i=1:ni
    
    % current number of needed frames 
    cnf = objectList(i).lft; 
    
    % start position of this object
    xystart = objectList(i).startpos;
    
    %generate particle's track
    if cnf > 1
        xyvecTraj = brownianMotion(2,diffCoef2D,cnf-1,1,1,confRad2D); %Brownian part
        if objectList(i).mType == 1 %linear part
            trackLin = brownianMotion(1,diffCoef1D,cnf-1,1);
            xyvecTraj(:,2) = xyvecTraj(:,2) + trackLin;
        end
    else
        xyvecTraj = zeros(1,2);
    end
    xyvecTraj = xyvecTraj + repmat(xystart,cnf,1); %add initial position
    
    % crop the trajectory (set those points that have wandered outside the
    % physical image to zero)
    xyCropTraj = cropTrajToImage(xyvecTraj,imSize);
    
    % start intensity of this object
    startint = objectList(i).startint;
    
    % current: assign random intensity from intVec, which implies that all
    % objects have the same intensity and that the distribution represents
    % frame-to-frame variation  
    intVecTraj = [ startint; intVec(1)+intVec(2)*randn(cnf-1,1) ];
    
    % however, don't allow negative intensities
    intVecTraj(intVecTraj<0) = 0;
    
    % and also set intensities to zero where x,y are zero
    intVecTraj( (xyCropTraj(:,1)==0) | (xyCropTraj(:,2)==0) ) = 0;
    
    % enter the positions and intensities into the MPM file
    sf = objectList(i).startframe; %startframe
    ef = sf + cnf - 1; %endframe
    simMPM(i,sf*3-2:3:ef*3-2) = xyCropTraj(:,1);
    simMPM(i,sf*3-1:3:ef*3-1) = xyCropTraj(:,2);
    simMPM(i,sf*3:3:ef*3) = intVecTraj;

end

    

end % of function



%% Subfunction 1

function xyCrop = cropTrajToImage(xy,imSize)

sl = size(xy,1);
xyCrop = xy;

for i=1:sl
    x0 = xy(i,1);
    y0 = xy(i,2);
    
    if (x0<1) || (y0<1) || (x0>imSize(1)) || (y0>imSize(2))
        xyCrop(i,:) = [0 0];
    end
    
end

end % of subfunction
