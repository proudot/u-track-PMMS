function [simMPM]=simulateMPMlftbeh(imsize, nump, lftDist, numf, intVec, motionVec);
% simulateMPMlftbeh simulates a distribution of objects with specified
% lifetimes and specified types of behavior; the results are written into
% an MPM file (containing x,y, and intensity), which can later be
% transformed into a simulated movies with the function makeAiryImageFromMPM
%
% INPUT: 	imsize      = image size vector [sx,sy]
%           nump        = number of points 
%                         Important: this refers to the average density of
%                         objects per image, there are variations depending
%                         on the lifetime distribution
%           lftDist     = lifetime distribution vector (normalized 
%                         probability). the vector is 1-dimensional, as the
%                         length position automatically corresponds to the 
%                         number of frames - if e.g. all objects should
%                         have the same lifetime 10 frames, then lftDist
%                         should have the form [0 0 0 0 0 0 0 0 0 1]           
%           numf        = number of frames
%           intVec      = intensity vector [average std] std refers to the
%                         pit-to-pit variation in intensity
%           motionVec   = motion vector, can be used to choose between
%                         different modulues for determining frame-to-frame
%                         motion or to set motion parameters


%% ========================================================================
%   intialize variables
%% ========================================================================

% maximum lifetime
pl = length(lftDist);
% x-length of image
lx = imsize(1);
% y-length of image
ly = imsize(2);
% when we assign start positions inside the image, we use a frame, so that
% the objects don't wander off the image too easily - this constraint can 
% be changed as necessary
frame = 5;


%% ========================================================================
%   determine number of iterations based on number of objects and frames
% =========================================================================

% expectancy value for lifetime
ex_lft = sum([1:pl].*shiftdim(lftDist)');

% the necessary number of simulated objects to reach the required specified
% density of objects per image is approximately
ni = round(nump*numf/ex_lft); 


%% ========================================================================
%   determine appropriate lifetimes for the speicifed number of objects
% =========================================================================

% loop over ni
i=1;
while i<=ni
    
    % randomly select object's starting frame, allow the search to go back
    % 100 frames before start of the movie
    startframe = round((100+numf)*rand(1))-100;
    
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
        vis_endframe   = min(numf,endframe);
        vis_lifetime   = vis_endframe - vis_startframe + 1;
        
        % enter values into object list
        objectList(i).startframe = vis_startframe;
        objectList(i).lft = vis_lifetime;
        
        % assign this object a random start position based on the intvec
        % distribution
        px = frame + 1 + ((lx-2*frame-1)*rand(1));
        py = frame + 1 + ((ly-2*frame-1)*rand(1));
        objectList(i).startpos = [px py];
        
        % assign this object a random start intensity
        objectList(i).startint = intVec(1)+intVec(2)*randn(1);

        % enter the actual lifetine into the allLFT vector, which allows us
        % to check the simulated lifetime distribution later on if
        % necessary
        allLFT(i) = plft;
        
        i = i+1;
    end
end

% CHECKPOINT: if desired, compare the resulting lifetime distribution to 
% the original probability density at this point
allLFThist = hist(allLFT,[1:pl]);
allLFTnorm = allLFThist/sum(allLFThist);
% comment/uncomment the following paragraph to display the simulated
% lifetime distribution
figure;
bar(allLFTnorm);
hold on;
plot(lftDist,'r.-');


%% ========================================================================
%   create trajectories for all objects in the list, and enter them into
%   the MPM files
% =========================================================================

% initialize mpm file 
simMPM = zeros(ni,3*numf);

for i=1:ni
    % start position of this object
    xystart = objectList(i).startpos;
    % current number of needed frames 
    cnf = objectList(i).lft; 
    
    % ====================================================================
    % create trajectory of the desired length for this object, using the
    % function and/or parameters specified in motionVec
    % ====================================================================

    [xyvecTraj] = createTraj(xystart,cnf,motionVec);
    
    % crop the trajectory (set those points that have wandered outside the
    % physical image to zero)
    [xyCropTraj] = cropTrajToImage(xyvecTraj,[imsize]);
    
    
    % ====================================================================
    % create intensity vector of the desired length for this object, using 
    % specified intVec (different kinds of intensity variation can be
    % implemented here)
    % ====================================================================
    
    % start intensity of this object
    startint = objectList(i).startint;
    % current: assign random intensity from intVec, which implies that all
    % objects have the same intensity and that the distribution represents
    % frame-to-frame variation  
    intVecTraj = [ startint; intVec(1)+intVec(2)*randn(cnf-1,1) ];
    % however, don't allow negative intensities
    intVecTraj(find(intVecTraj<0)) = 0;
    % and also set intensities to zero where x,y are zero
    zeroPos = find( (xyCropTraj(:,1)==0) | (xyCropTraj(:,2)==0) );
    intVecTraj(zeroPos) = 0;
    
    
    % ====================================================================
    % enter the positions and intensities into the MPM file
    % ====================================================================
    
    % startframe
    sf = objectList(i).startframe;
    % endframe
    ef = sf + cnf - 1;
    
    simMPM(i,sf*3-2:3:ef*3-2) = xyCropTraj(:,1);
    simMPM(i,sf*3-1:3:ef*3-1) = xyCropTraj(:,2);
    simMPM(i,sf*3:3:ef*3) = intVecTraj;
    
    
end

    

end % of function




%% ========================================================================
%
%       SUBFUNCTIONS
%
% =========================================================================

%% ========================================================================
%   implement whatever kind of motion you want here, using e.g. Gunnar's
%   code to generate trajectories composed of different motion types
% =========================================================================

function [xyTraj] = createTraj(xystart,nf, motionVec);

xyTraj = zeros(nf,2);
xyTraj(1,:) = xystart;

xycurr = xyTraj(1,:);
% current version: simple random walk
diff = motionVec;
for k=2:nf
    dx = diff*randn(1);
    dy = diff*randn(1);
    xynew = xycurr + [dx dy];
    xycurr = xynew;
    xyTraj(k,:) = xynew;
end

end % of subfunction





function [xyCrop]=cropTrajToImage(xy,imsize);

[sl,s2] = size(xy);
xyCrop = xy;

for i=1:sl
    x0 = xy(i,1);
    y0 = xy(i,2);
    
    if (x0<1) | (y0<1) | (x0>imsize(1)) | (y0>imsize(2))
        xyCrop(i,:) = [0 0];
    end
    
end

end % of subfunction
