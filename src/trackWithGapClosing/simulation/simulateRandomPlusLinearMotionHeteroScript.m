imSize=[200 200];
density=3;
numP=density*(imSize(1)*imSize(2))/400
%numP=50;
numF=100;
sigmaLft=10;
meanLft=50;
x =1:numF;
lftDist=exp(-(x-meanLft).^2/(2*sigmaLft^2))/sqrt(2*pi*sigmaLft^2);
intVec=[500 30];
bitDepth=8;

clear motionParam
motionParam.diffCoef2D=[1 3];
motionParam.confRad2D=[0.5 4];
motionParam.speed1D=[1 6];
motionParam.probTypeSwitch=0.05;
motionParam.probVelSwitch=0.2;
motionParam.fracType=[0 0 0 0 0 0 1]

motionParam.probMS=[];

[simMPM,tracksSim] = simulateRandomPlusLinearMotionHetero(imSize,numP,lftDist,numF,intVec,motionParam);
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

%% data translation
% detection GT to movie info data structure
% MovieInfo data structure allows to assign an implicite ID to each spot
% We have to store this tracks in the ID form to link kalman filter info to real GT tracks property.

movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numF,1);
for i = 1: numF
    xCoord=simMPM(simMPM(:,3*i-2)~=1.,3*i-2);
    movieInfo(i).xCoord=[ xCoord 0.05*ones(size(xCoord))];
    yCoord=simMPM(simMPM(:,3*i-1)~=1.,3*i-1);
    movieInfo(i).yCoord=[ yCoord 0.05*ones(size(yCoord))];    
    amp=simMPM(simMPM(:,3*i)~=0.,3*i)/(2^bitDepth -1);
    movieInfo(i).amp=[ amp 0.0003*ones(size(amp))];        
end

bgav=100;
bgnoise=30;
sigma=1.;
rad=4;

saveVar=1;
saveFolder='/home/pr120/code/heterogenous_track/data_set/vimentin/simulation/sandbox/';
imsize=imSize;
% solving strange (x,y) reversing issue with makeAiryImageFromMPM
revSimMPM=simMPM;                       
revSimMPM(:,1:3:end)=simMPM(:,2:3:end);
revSimMPM(:,2:3:end)=simMPM(:,1:3:end);
[imageStack]=makeAiryImageFromMPM(revSimMPM,bgav,bgnoise,sigma,imsize,rad,saveVar,saveFolder);

imageStacku8=uint8(imageStack);
implay(imageStacku8);

% makeAiryImageFromMPM makes an image from the point distribution specified
% in the mpm file, using specified values for amplitudes and noise
%
% SYNOPSIS:
% [imageStack]=makeAiryImageFromMPM(mpm,amps,bg,bgnoise,sigma,imsize,rad)
%
% INPUT     :   trackInfo = MPM-like file with location and intensities of 
%                           objects in successive columns:
%                           x1,y1,i1,x2,y2,i2,...
%               bgav      = average background level (can be either single 
%                           value or, for background inequality, an image 
%                           of the same size as imsize). In counts, e.g.
%                           assuming a 16-bit camera.
%               bgnoise   = std of backgound noise. In counts, e.g.
%                           assuming a 16-bit camera.
%               sigma     = width of point-spread function (in pixels),
%                           where sigma = (1/3)*(radius of Airy disc)
%               imsize    = size of image [sx,sy]
%               rad       = radius used for the generation of the Airy disc
%                           for each object, in increments of sigma 
%                           ( should ideally be >=3 )
%               saveVar   = variable that indicates whether or not tif files
%                           are saved to file (1/0)
%               saveFolder  = (optional) folder name
%
% OUTPUT    :   imageStack
%
% REMARKS
%



movieParam.imageDir = saveFolder;
movieParam.filenameBase = 'simulImage'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum =numF; %number of last image in movie
movieParam.digits4Enum = 3; %number of digits used for frame enumeration (1-4).

saveResults.dir = [ movieParam.imageDir 'res/' ];

revImageDir = [movieParam.imageDir 'rev/'];