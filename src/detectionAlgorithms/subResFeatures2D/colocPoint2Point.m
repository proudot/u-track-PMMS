function [fracPointsOverlap,pValue,colocType,numPoints12,fracPointsOverlapRand] = ...
    colocPoint2Point(pointImageRef,pointImageColoc,pointDetectInput,...
    colocRadius,numRepTestRand,doPlot)
%COLOCPOINT2POINT estimates the colocalization fraction and its significance between two sets of points
%
%SYNOPSIS [fracPointsOverlap,pValue,colocType,numPoints12,fracPointsOverlapRand] = ...
%    colocPoint2Point(pointImageRef,pointImageColoc,pointDetectInput,...
%    colocRadius,numRepTestRand,doPlot)
%
%INPUT  pointImageRef   : Image of reference objects to be colocalized
%                         against.
%       pointImageColoc : Image of objects to be colocalized with
%                         reference objects.
%       pointDetectInput:
%                       EITHER
%                         2x1 array containing the alpha-value for local
%                         maxima detection in the two images.
%                       OR
%                         2x1 structure array containing the fields
%                         .xCoord and .yCoord for both images 
%                         (i.e. concatenated movieInfo, the output of
%                         detectSubResFeatures2D_StandAlone).
%                       Optional. Default: alphaLocMax = 0.05 for both
%                       images.
%       colocRadius     : Radius to look for colocalization, in pixels.
%                         Optional. Default: 2.
%       numRepTestRand  : Number of times to sample a random distribution
%                         of point positions (with the same number of
%                         points) to generate a distribution of perchance
%                         colocalization fractions and test the
%                         significance of the observed colocalization
%                         fraction.
%                         Optional. Default: 500.
%       doPlot          : 1 to plot results, 0 otherwise. If 1, the code
%                         will make 3 plots - one showing the detected
%                         points on top of pointImageRef, one showing the
%                         detected lines on top of pointImageColoc, and one
%                         showing the detected points on top of each other.
%                         In the thrid figure, points considered
%                         colocalized will be in green, those not
%                         colocalized will be in red.
%                         Optional. Default: 1.
%
%OUTPUT fracPointsOverlap    : Fraction of points in pointImageColoc
%                              colocalizing with points in pointImageRef.
%                              Fraction is relative to total number of
%                              points in pointImageColoc.
%       pValue               : P-value indicating the significance of the
%                              observed colocalization fraction. P-value
%                              calculated from the distribution of
%                              perchance colocalization fractions.
%       colocType            : 'attract'/'repel' to indicate that fraction
%                              of colocalizatoin indicates attraction/repulsion.
%       numPoints12          : Number of detected points in both images.
%       fracPointsOverlapRand: Distribution of perchance colocalization
%                              fractions.
%
%Khuloud Jaqaman, April 2011

%% Input

if nargin < 2
    disp('--colocPoint2Point: Please input images')
    return
end

%rename images for convenience
pointImage1 = pointImageRef;
pointImage2 = pointImageColoc;

%point detection parameters
if nargin < 3 || isempty(pointDetectInput)
    alphaLocMax1 = 0.05;
    alphaLocMax2 = 0.05;
else
    if isstruct(pointDetectInput)
        movieInfo1 = pointDetectInput(1);
        movieInfo2 = pointDetectInput(2);
        alphaLocMax1 = [];
        alphaLocMax2 = [];
    else
        alphaLocMax1 = pointDetectInput(1);
        alphaLocMax2 = pointDetectInput(2);
    end
end

%colocalization radius
if nargin < 4 || isempty(colocRadius)
    colocRadius = 2;
else
    colocRadius = ceil(colocRadius);
end

%number of repetitions to estimate significance
if nargin < 5 || isempty(numRepTestRand)
    numRepTestRand = 500;
end

%plot or not
if nargin < 6 || isempty(doPlot)
    doPlot = 1;
end

%get image size
[imgSizeX,imgSizeY] = size(pointImageRef);

%% Get object coordinates

for iImage = 1 : 2
    
    %read this image's information
    eval(['pointImage = pointImage' num2str(iImage) ';'])
    eval(['alphaLocMax = alphaLocMax' num2str(iImage) ';'])
    
    if ~isempty(alphaLocMax) %detect the objects
        
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
        
    else %or extract positions from movieInfo
        
        eval(['movieInfo = movieInfo' num2str(iImage) ';'])
        positions = round([movieInfo.yCoord(:,1) movieInfo.xCoord(:,1)]);
        
    end
    
    %convert positions from a 2D index to a 1D index
    positions1D = (positions(:,2)-1)*imgSizeX + positions(:,1);
    
    %get total number of detected points
    numPoints12(iImage) = size(positions,1);
    
    %store this image's detection results
    eval(['positions' num2str(iImage) ' = positions;'])
    eval(['positions' num2str(iImage) '_1D = positions1D;'])
    
end

%% Colocalization

%place circles around objects in reference image
image1Circles = zeros(imgSizeX,imgSizeY);
image1Circles(positions1_1D) = 1;
SE = strel('disk',colocRadius);
image1Circles = imdilate(image1Circles,SE);

%look at fraction of points in image 2 colocalizing with points in image 1
indxColoc = find(image1Circles(positions2_1D)==1);
indxNotColoc = setdiff(1:numPoints12(2),indxColoc);
fracPointsOverlap = length(find(image1Circles(positions2_1D)==1))/numPoints12(2);
fracPointsOverlap = length(indxColoc)/numPoints12(2);

%% Test significance of colocalization by randomizing positions

%get how many times to repeat calculation
numRep = numRepTestRand;

%generate numPoints with random positions, numRep times
randomPos = ceil(rand(numPoints12(2),numRep)*imgSizeX*imgSizeY);

%get number of randomly distributed points that colocalize with objects in
%reference image
fracPointsOverlapRand = NaN(numRep,1);
for iRep = 1 : numRep
    fracPointsOverlapRand(iRep) = length(find(image1Circles(randomPos...
        (:,iRep))==1))/numPoints12(2);
end

%calculate mean random colocalization fraction
medianColocRand = median(fracPointsOverlapRand);

%calculate p-value and type of colocalization
if fracPointsOverlap >= medianColocRand
    pValue = length(find(fracPointsOverlapRand>=fracPointsOverlap))/numRep;
    colocType = 'attract';
else
    pValue = length(find(fracPointsOverlapRand<=fracPointsOverlap))/numRep;
    colocType = 'repel';
end
pValue = min(pValue,0.5);
    
%% Plot if requested

if doPlot

    %plot detected points on top of pointImageRef
    figure('Name','Reference objects')
    imshow(pointImageRef,[])
    hold on
    plot(positions1(:,2),positions1(:,1),'r.')
    
    %plot detected points on top of pointImageColoc
    figure('Name','Objects to colocalize')
    imshow(pointImageColoc,[])
    hold on
    plot(positions2(:,2),positions2(:,1),'r.')
    
    %plot detected points on top of each other
    figure('Name','Colocalization')
    imshow(image1Circles,[])
    hold on
    plot(positions2(indxColoc,2),positions2(indxColoc,1),'g.')
    plot(positions2(indxNotColoc,2),positions2(indxNotColoc,1),'r.')
    
end

%% ~~~ the end ~~~


