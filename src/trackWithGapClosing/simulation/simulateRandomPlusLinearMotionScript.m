imSize=[150 150];
numP=100;
numF=50;
sigmaLft=10;
meanLft=20;
x =1:numF;
lftDist=exp(-(x-meanLft).^2/(2*sigmaLft^2))/sqrt(2*pi*sigmaLft^2);
intVec=[500 30];
motionParam.diffCoef2D=[0.1 1];
motionParam.confRad2D=[0.5 2];
motionParam.speed1D=[4 10];
motionParam.probVelSwitch=0.8;
%motionParam.fracType=[0.1 0.7 0.2];
motionParam.fracType=[0.5 0.4 0.1];

motionParam.probMS=[];

[simMPM,tracksSim] = simulateRandomPlusLinearMotion(imSize,numP,lftDist,numF,intVec,motionParam);
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
simMPM

bgav=100;
bgnoise=30;
sigma=1.;
rad=4;

saveVar=1;
saveFolder='/home/pr120/code/heterogenous_track/data_set/vimentin/simulation/sandbox/';
imsize=imSize;
[imageStack]=makeAiryImageFromMPM(simMPM,bgav,bgnoise,sigma,imsize,rad,saveVar,saveFolder);
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
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% DATE: 28-Feb-2006
% last modified
% DATE: 05-Oct-2007
%
%

