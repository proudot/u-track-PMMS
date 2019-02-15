function [movieInfo,meanBg,tracksStatic] = rawIntFromMovieInfo(movieInfo,psfSigma,imageFilePath,static)
%RAWINTFROMMOVIEINFO reads raw image intensities at locations specified in movieInfo
%
%SYNOPSIS movieInfo = rawIntFromMovieInfo(movieInfo,psfSigma,imageFilePath)
%
%INPUT  
%       movieInfo    : Array of size equal to the number of frames
%                      in movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%             .yCoord      : Image coordinate system y-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D. Default: zeros.
%             .zCoord      : Image coordinate system z-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D or 2D. Default: zeros.
%             .amp         : Amplitudes of PSFs fitting detected features. 
%                            1st column for values and 2nd column 
%                            for standard deviations.
%       psfSigma     : Point Spread Function sigma in pixels.
%                      Optional. Default: 1 pixel. 
%       imageFilePath: File name, including path, of tiff file storing
%                      images. Code assumes that images are stored in a
%                      multi-image tiff file.
%                      Optional. If not input, user will be prompted to
%                      specifiy file.
%       static       : Flag with value 0 if object locations are not
%                      coupled between frames, or value 1 to take
%                      the locations in frame 1 and read those same
%                      locations in all subsequent frames.
%                      Optional. Default: 0.
%
%OUTPUT movieInfo    : Same as input, but with added fields "intRaw" and
%                      "intRawMinusBg." Note that the reported intensity is
%                      the average intensity (with or without background
%                      subtraction) over an area of
%                      5*psfSigsm x 5*psfSigma (so for psfSigma=1, a 5x5 area).
%       meanBg       : Mean background value per frame.
%
%Khuloud Jaqaman, July 2014

%% Input

if nargin < 2 || isempty(psfSigma)
    psfSigma = 1;
end

if nargin < 3 || isempty(imageFilePath)
    [fName,dirName] = uigetfile('*.tif','Please specify image file');
    imageFilePath = fullfile(dirName,fName);
end

if nargin < 4 || isempty(static)
    static = 0;
end

%find number of frames in movie
numFrames = length(movieInfo);

%find image size
image = double(imread(imageFilePath,1));
imSize = size(image);

%% Read intensity

%modify movieInfo if static = 1
if static
    for iFrame = 2 : numFrames
        movieInfo(iFrame) = movieInfo(1);
        movieInfo(iFrame).amp = [];
    end
end

%determine patch size for reading raw intensity
patchLength = floor(2.5*psfSigma);

%go over all frames and read raw intensity
meanBg = zeros(numFrames,1);
for iFrame = 1 : numFrames

    %read in image
    image = double(imread(imageFilePath,iFrame));
    
    %get all positions
    xCoord = round(movieInfo(iFrame).xCoord);
    yCoord = round(movieInfo(iFrame).yCoord);
    
    %get number of objects
    numObj = size(xCoord,1);
    intRaw = zeros(numObj,1);
    
    %make default backgroun image
    imageBg = image;
    
    %if there are objects in this frame
    if numObj > 0
        
        xCoord = xCoord(:,1);
        yCoord = yCoord(:,1);
        
        %go over all positions and read intensity
        for iObj = 1 : numObj
            minX = max(1,xCoord(iObj)-patchLength);
            minY = max(1,yCoord(iObj)-patchLength);
            maxX = min(imSize(1),xCoord(iObj)+patchLength);
            maxY = min(imSize(2),yCoord(iObj)+patchLength);
            imgObj = image(minY:maxY,minX:maxX);
            intRaw(iObj) = mean(imgObj(:));
        end
        
        %get only background areas
        linIndx = sub2ind(imSize,yCoord,xCoord);
        bgMask = zeros(imSize);
        bgMask(linIndx) = 1;
        bgMask = imdilate(bgMask,strel('square',patchLength*2+1));
        bgMask = 1 - bgMask;
        imageBg = image .* bgMask;
        
    end
    
    meanBg(iFrame) = mean(imageBg(:));
    
    %substract background intensity from raw intensity
    intRawMinusBg = intRaw - meanBg(iFrame);
    
    %store results for output
    movieInfo(iFrame).intRaw = intRaw;
    movieInfo(iFrame).intRawMinusBg = intRawMinusBg;
    
end

%% Tracks in case of static

if static
    
    %track matrix
    tracksMat = repmat([xCoord yCoord zeros(numObj,6)],1,numFrames);
    for iFrame = 1 : numFrames
        tracksMat(:,8*(iFrame-1)+4) = movieInfo(iFrame).intRaw;
    end
    
    %convert to structure format
    tracksStatic = repmat(struct('tracksCoordAmpCG',[],'tracksFeatIndxCG',[],'seqOfEvents',[]),numObj,1);
    for iObj = 1 : numObj
        tracksStatic(iObj).tracksCoordAmpCG = tracksMat(iObj,:);
        tracksStatic(iObj).tracksFeatIndxCG = iObj*ones(1,numFrames);
        tracksStatic(iObj).seqOfEvents = [1 1 1 NaN; numFrames 2 1 NaN];
    end
    
else
    
    tracksStatic = [];
    
end

%%%%% ~~ the end ~~ %%%%%
