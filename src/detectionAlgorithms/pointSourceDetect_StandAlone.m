function [movieInfo] =  pointSourceDetect_StandAlone(movieParam,saveResults,verbose,sigma,bitDepth,varargin)
%pointSourceDetect_StandAlone is a Q&D pointSourceDetect wrapper for U-track.
%
%SYNOPSIS [movieInfo,exceptions,localMaxima,background,psfSigma] = ...
%    pointSourceDetect_StandAlone(movieParam,detectionParam,saveResults)
%
%INPUT  movieParam    : Structure with fields
%           .imageDir     : Directory where images are stored.
%           .filenameBase : Filename base.
%           .firstImageNum: Numerical index of first image in movie.
%           .lastImageNum : Numerical index of last image in movie.
%           .digits4Enum  : Number of digits used to enumerate frames.


%% Output

movieInfo = [];
pstruct=[];

%% input + init

%get movie parameters
hasImageDir = isfield(movieParam,'imageDir');
if hasImageDir
    imageDir = movieParam.imageDir;
    filenameBase = movieParam.filenameBase;
    digits4Enum = movieParam.digits4Enum;
else
    channel = movieParam.channel;
end
firstImageNum = movieParam.firstImageNum;
lastImageNum = movieParam.lastImageNum;

if hasImageDir
    %store the string version of the numerical index of each image
    enumString = repmat('0',lastImageNum,digits4Enum);
    formatString = ['%0' num2str(digits4Enum) 'i'];
    for i=1:lastImageNum
        enumString(i,:) = num2str(i,formatString);
    end
end


%% General image information

%get image indices and number of images
imageIndx = firstImageNum : lastImageNum;
numImagesRaw = lastImageNum - firstImageNum + 1; %raw images

if hasImageDir
    %read first image and get image size
    if exist([imageDir filenameBase enumString(imageIndx(1),:) '.tif'],'file')
        imageTmp = imread([imageDir filenameBase enumString(imageIndx(1),:) '.tif']);
    else
        disp('First image does not exist! Exiting ...');
        return
    end
    [imageSizeX,imageSizeY] = size(imageTmp);
    clear imageTmp
else
    imageSizeX=channel.owner_.imSize_(1);
    imageSizeY=channel.owner_.imSize_(2);
end

%check which images exist and which don't
imageExists = true(numImagesRaw,1);
if hasImageDir
    for iImage = 1 : numImagesRaw
        if ~exist([imageDir filenameBase enumString(imageIndx(iImage),:) '.tif'],'file')
            imageExists(iImage) = 0;
        end
    end
end

%initialize movieInfo
clear movieInfo
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numImagesRaw,1);

for iImage = 1 : numImagesRaw
    if(verbose)
       disp(['Processing frame ' int2str(iImage)] ) ;
    end
    currentImage=double(imread([imageDir filenameBase ...
                        enumString(imageIndx(iImage),:) '.tif']));

    currentImage = currentImage / (2^bitDepth-1); % normalize
                                
    [pstruct]=pointSourceDetection(currentImage,sigma,varargin{:});
    if(~(isempty(pstruct)))
        featuresInfo.xCoord=[ pstruct.x' pstruct.x_pstd'];
        featuresInfo.yCoord=[ pstruct.y' pstruct.y_pstd'];
        featuresInfo.amp=[ pstruct.A' pstruct.A_pstd'];
    else 
        featuresInfo.xCoord=[ ];
        featuresInfo.yCoord=[ ];
        featuresInfo.amp=[ ];
    end
    movieInfo(iImage) = featuresInfo;
end    


%save results
if isstruct(saveResults)
    save([saveResults.dir filesep saveResults.filename],'movieParam', 'movieInfo');
end
