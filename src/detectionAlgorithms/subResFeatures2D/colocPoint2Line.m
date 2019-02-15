function [fracPointsOnLines,pValue,numPoints,fracPointsOnLinesRand] = ...
    colocPoint2Line(lineImage,pointImage,lineDetectInput,pointDetectInput,...
    numRepTestRand,doPlot)
%COLOCPOINT2LINE estimates the colocalization fraction and its significance between a set of points and a set of lines
%
%SYNOPSIS [fracPointsOnLines,pValue,numPoints,fracPointsOnLinesRand] = ...
%    colocPoint2Line(lineImage,pointImage,lineDetectInput,pointDetectInput,...
%    numRepTestRand,doPlot)
%
%INPUT  lineImage       : Image of line features.
%       pointImage      : Image of point features.
%       lineDetectInput : Structure of parameters needed for imLineDetect
%                         that detects the lines in lineImage.
%                         Contains the fields:
%           .scales         : Same as input variable "scales" in
%                             imLineDetect. See imLineDetect for help.
%                             Optional. Default: [0.5:0.5:4].
%           .lineType       : Same as input variable "lineType" in
%                             imLineDetect. See imLineDetect for help.
%                             Optional. Default: 1.
%           .conf           : Same as input variable "conf" in
%                             imLineDetect. See imLineDetect for help.
%                             Optional. Default: 0.9.
%           .dilateW        : Width of square to be used to dilate the line
%                             response in lineImage for the final
%                             construction of the detected lines.
%                             Optional. Default: 1 (i.e. no dilation).
%                         Whole structure optional.
%       pointDetectInput: Structure of parameters needed for detecting the
%                         local maxima representing the point features in
%                         pointImage.
%                         Contains the fields:
%           .alphaLocMax    : Alpha-value for testing the significance of
%                             local maxima against background.
%                             Optional. Default: 0.05.
%                         Whole structue optional.
%       numRepTestRand  : Number of times to sample a random distribution
%                         of point positions (with the same number of
%                         points) to generate a distribution of perchance
%                         colocalization fractions and test the
%                         significance of the observed colocalization
%                         fraction.
%                         Optional. Default: 200.
%       doPlot          : 1 to plot results, 0 otherwise. If 1, the code
%                         will make 3 plots - one showing the detected
%                         points on top of pointImage, one showing the
%                         detected lines on top of lineImage, and one
%                         showing the detected points on top of the
%                         detected lines.
%                         Optional. Default: 1.
%
%OUTPUT fracPointsOnLines    : Fraction of point objects colocalizing with
%                              line objects.
%       pValue               : P-value indicating the significance of the
%                              observed colocalization fraction. P-value
%                              estimated from the distribution of
%                              perchance colocalization fractions.
%       numPoints            : Number of detected point objects.
%       fracPointsOnLinesRand: Distribution of perchance colocalization
%                              fractions.
%
%Khuloud Jaqaman, Summer 2008

%% input

if nargin < 2
    disp('--colocPoint2Line: Please input at least line object image and point object image')
    return
end

%line detection parameters
if nargin < 3 || isempty(lineDetectInput)
    scales = 0.5:0.5:4;
    lineType = 1;
    conf = 0.9;
    dilateW = 1;
else
    if isfield(lineDetectInput,'scales')
        scales = lineDetectInput.scales;
    else
        scales = 0.5:0.5:4;
    end
    if isfield(lineDetectInput,'lineType')
        lineType = lineDetectInput.lineType;
    else
        lineType = 1;
    end
    if isfield(lineDetectInput,'conf')
        conf = lineDetectInput.conf;
    else
        conf = 0.9;
    end
    if isfield(lineDetectInput,'dilateW')
        dilateW = lineDetectInput.dilateW;
    else
        dilateW = 1;
    end
end

%point detection parameters
if nargin < 4 || isempty(pointDetectInput)
    alphaLocMax = 0.05;
else
    if isfield(pointDetectInput,'alphaLocMax')
        alphaLocMax = pointDetectInput.alphaLocMax;
    else
        alphaLocMax = 0.05;
    end
end

%number of repetitions to estimate significance
if nargin < 5 || isempty(numRepTestRand)
    numRepTestRand = 200;
end

%plot or not
if nargin < 6 || isempty(doPlot)
    doPlot = 1;
end

%get image size
[imgSizeX,imgSizeY] = size(pointImage);

%% line detection

%detect lines in lineImage
img.data = lineImage;
img.perm = 'M';
respRaw = imLineDetect(img,scales,lineType,conf);

%keep only significant lines - HOW???
% respVal = respRaw(respRaw~=0);
% respThresh = prctile(respVal,5);
respThresh = 0;
resp = respRaw;
resp(resp<=respThresh) = 0;
resp(resp>respThresh) = 1;

%dilate resp by requested width
SE = strel('square',dilateW);
respDilate = imdilate(resp,SE);

%% point detection

%filter pointImage
pointImageF = filterGauss2D(pointImage,1);

%get background information
[bgMean,bgStd] = spatialMovAveBG(pointImageF,size(pointImageF,1),size(pointImageF,2));

%find local maxima in filtered image
fImg = locmax2d(pointImageF,[3 3],1);

%get positions and amplitudes of local maxima
[localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);
localMax1DIndx = find(fImg(:));

%get background values corresponding to local maxima
bgMeanMax = bgMean(localMax1DIndx);
bgStdMax = bgStd(localMax1DIndx);

%calculate the p-value corresponding to the local maxima's amplitudes
%assume that background intensity in filtered image is normally
%distributed with mean bgMeanMax and standard deviation bgStdMax
pValue = 1 - normcdf(localMaxAmp,bgMeanMax,bgStdMax);

%retain only those maxima with significant amplitude
keepMax = find(pValue < alphaLocMax);
localMaxPosX = localMaxPosX(keepMax);
localMaxPosY = localMaxPosY(keepMax);
positions = [localMaxPosX localMaxPosY];

%convert positions from a 2D index to a 1D index
positions1D = (positions(:,2)-1)*imgSizeX + positions(:,1);

%get total number of detected points
numPoints = size(positions,1);

%% colocalization

%get number of points that colocalize with lines
fracPointsOnLines = length(find(respDilate(positions1D)==1))/numPoints;

%% test colocalization with a random distribution

%get how many times to repeat calculation
numRep = numRepTestRand;

%generate numPoints with random positions, numRep times
randomPos = ceil(rand(numPoints,numRep)*imgSizeX*imgSizeY);

%get number of randomly distributed points that colocalize with lines
fracPointsOnLinesRand = NaN(numRep,1);
for iRep = 1 : numRep
    fracPointsOnLinesRand(iRep) = length(find(respDilate(randomPos...
        (:,iRep))==1))/numPoints;
end

%% assuming a normal distribution, get observed value's p-value

meanRand = mean(fracPointsOnLinesRand);
stdRand = std(fracPointsOnLinesRand);

if fracPointsOnLines > meanRand
    pValue = 1 - normcdf(fracPointsOnLines,meanRand,stdRand);
else
    pValue = normcdf(fracPointsOnLines,meanRand,stdRand);
end

%% plot if requested
if doPlot

    %plot detected points on top of pointImage
    plotImageWithFeatures(pointImage,[positions(:,2) positions(:,1)]);
    
    %plot detected lines on top of lineImage
    [linePosX,linePosY] = find(resp~=0);
    plotImageWithFeatures(lineImage,[linePosY linePosX]);
    
    %plot detected points on top of detected lines
    plotImageWithFeatures(respDilate,[positions(:,2) positions(:,1)]);
    
end

%% ~~~ the end ~~~


