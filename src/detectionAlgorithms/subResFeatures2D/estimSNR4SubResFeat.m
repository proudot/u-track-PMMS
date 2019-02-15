function signal2noiseRatio = estimSNR4SubResFeat(movieInfo,movieParam,bitDepth)
%ESTIMSNR4SUBRESFEAT estimates the SNR in movies with a (possibly) spatially heterogeneous background
%
%SYNOPSIS signal2noiseRatio = estimSNR4SubResFeat(movieInfo,movieParam,bitDepth)
%
%INPUT  movieInfo     : Structure array of length = number of frames in
%                       movie, containing the fields:
%           .xCoord       : Image coordinate system x-coordinate of detected
%                           features [x dx] (in pixels).
%           .yCoord       : Image coordinate system y-coordinate of detected
%                           features [y dy] (in pixels).
%           .amp          : Amplitudes of PSFs fitting detected features [a da].
%       movieParam    : Structure with fields
%           .imageDir     : Directory where images are stored.
%           .filenameBase : Filename base.
%           .firstImageNum: Numerical index of first image in movie.
%           .lastImageNum : Numerical index of last image in movie.
%           .digits4Enum  : Number of digits used to enumerate frames.
%       bitDepth      : Camera bit depth.
%
%OUTPUT signal2noiseRatio: Number of features - by - number of frames
%                          array showing signal to noise ratio of all
%                          features in all frames (SNR = signal amplitude
%                          above background / local background std).
%
%Khuloud Jaqaman, January 2008

%% input

%get movie parameters
imageDir = movieParam.imageDir;
filenameBase = movieParam.filenameBase;
firstImageNum = movieParam.firstImageNum;
lastImageNum = movieParam.lastImageNum;
digits4Enum = movieParam.digits4Enum;

%store the string version of the numerical index of each image
enumString = getStringIndx(digits4Enum);

%get number of features in each image
numFeat = zeros(length(movieInfo),1);
for iImage = 1 : length(movieInfo)
    numFeat(iImage) = size(movieInfo(iImage).xCoord,1);
end
numFeatMax = max(numFeat);

%% read images

%get image related parameters
imageIndx = firstImageNum : lastImageNum; %image indices
imageTmp = imread([imageDir filenameBase enumString(imageIndx(1),:) '.tif']); %first image
[imageSizeX,imageSizeY] = size(imageTmp); %image size
numImagesRaw = lastImageNum - firstImageNum + 1; %number of images
clear imageTmp

%read images and store them in array
imageRaw = zeros(imageSizeX,imageSizeY,numImagesRaw);
for iImage = 1 : numImagesRaw
    imageRaw(:,:,iImage) = imread([imageDir filenameBase enumString(imageIndx(iImage),:) '.tif']);    
end

%replace zeros with NaNs
%zeros result from cropping that leads to curved boundaries
imageRaw(imageRaw==0) = NaN;

%normalize images
imageRaw = double(imageRaw) / (2^bitDepth-1);

%% estimate background signal standard deviation

%initialize
bgStd = zeros(imageSizeX,imageSizeY,numImagesRaw);

%go over all images and get background signal statistics
for iImage = 1 : numImagesRaw
    [dummy,bgStd(:,:,iImage)] = spatialMovAveBG(imageRaw(:,:,iImage),...
        imageSizeX,imageSizeY);
end

%% divide feature amplitudes by local background signal std

signal2noiseRatio = NaN(numFeatMax,numImagesRaw);

%go over all images
for iImage = 1 : numImagesRaw

    %get feature coordinates, amplitudes and number of features in image
    xCoord = movieInfo(imageIndx(iImage)).xCoord(:,1);
    yCoord = movieInfo(imageIndx(iImage)).yCoord(:,1);
    amp = movieInfo(imageIndx(iImage)).amp(:,1);

    %get local background signal std
    bgStdLocal = bgStd(round(yCoord),round(xCoord),iImage);
    bgStdLocal = diag(bgStdLocal);
    
    %divide amplitudes by local background std
    snrValues = amp./bgStdLocal;
    
    %store in output matrix
    signal2noiseRatio(1:numFeat(iImage),iImage) = snrValues;

end


%% ~~~ the end ~~~

%% Subfunction 1

function enumString = getStringIndx(digits4Enum)

switch digits4Enum
    case 4
        enumString = repmat('0',9999,4);
        for i = 1 : 9
            enumString(i,:) = ['000' num2str(i)];
        end
        for i = 10 : 99
            enumString(i,:) = ['00' num2str(i)];
        end
        for i = 100 : 999
            enumString(i,:) = ['0' num2str(i)];
        end
        for i = 1000 : 9999
            enumString(i,:) = num2str(i);
        end
    case 3
        enumString = repmat('0',999,3);
        for i = 1 : 9
            enumString(i,:) = ['00' num2str(i)];
        end
        for i = 10 : 99
            enumString(i,:) = ['0' num2str(i)];
        end
        for i = 100 : 999
            enumString(i,:) = num2str(i);
        end
    case 2
        enumString = repmat('0',99,2);
        for i = 1 : 9
            enumString(i,:) = ['0' num2str(i)];
        end
        for i = 10 : 99
            enumString(i,:) = num2str(i);
        end
    case 1
        enumString = repmat('0',9,1);
        for i = 1 : 9
            enumString(i,:) = num2str(i);
        end
end

%% Subfunction 2

function [bgMean,bgStd] = spatialMovAveBG(imageLast5,imageSizeX,imageSizeY)

%the function in its current form assigns blocks of 11x11 pixels the
%same background values, for the sake of speed

%define pixel limits where moving average can be calculated
startPixelX = 16;
endPixelX = imageSizeX - 15;
startPixelY = 16;
endPixelY = imageSizeY - 15;

%allocate memory for output
bgMean = NaN(imageSizeX,imageSizeY);
bgStd = bgMean;

%go over all pixels within limits
for iPixelX = startPixelX : 11 : endPixelX
    for iPixelY = startPixelY : 11 : endPixelY
        
        %get local image
        imageLocal = imageLast5(iPixelX-15:iPixelX+15,iPixelY-15:iPixelY+15,:);
        
        %estimate robust mean and std
        %first remove NaNs representing cropped regions
        imageLocal = imageLocal(~isnan(imageLocal));
        if ~isempty(imageLocal)
            [bgMean1,bgStd1] = robustMean(imageLocal(:));
        else
            bgMean1 = NaN;
            bgStd1 = NaN;
        end
        
        %put values in matrix representing image
        bgMean(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5) = bgMean1;
        bgStd(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5) = bgStd1;
        
    end
end

%find limits of actual pixels filled up above
firstFullX = find(~isnan(bgMean(:,startPixelY)),1,'first');
lastFullX = find(~isnan(bgMean(:,startPixelY)),1,'last');
firstFullY = find(~isnan(bgMean(startPixelX,:)),1,'first');
lastFullY = find(~isnan(bgMean(startPixelX,:)),1,'last');

%patch the rest
for iPixelY = firstFullY : lastFullY
    bgMean(1:firstFullX-1,iPixelY) = bgMean(firstFullX,iPixelY);
    bgMean(lastFullX+1:end,iPixelY) = bgMean(lastFullX,iPixelY);
    bgStd(1:firstFullX-1,iPixelY) = bgStd(firstFullX,iPixelY);
    bgStd(lastFullX+1:end,iPixelY) = bgStd(lastFullX,iPixelY);
end
for iPixelX = 1 : imageSizeX
    bgMean(iPixelX,1:firstFullY-1) = bgMean(iPixelX,firstFullY);
    bgMean(iPixelX,lastFullY+1:end) = bgMean(iPixelX,lastFullY);
    bgStd(iPixelX,1:firstFullY-1) = bgStd(iPixelX,firstFullY);
    bgStd(iPixelX,lastFullY+1:end) = bgStd(iPixelX,lastFullY);
end

