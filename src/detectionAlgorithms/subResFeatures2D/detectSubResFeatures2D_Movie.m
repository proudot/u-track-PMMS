function [movieInfo,emptyFrames,framesFailed,errFlag] = detectSubResFeatures2D_Movie(...
    movieParam,detectionParam,saveResults)
%DETECTSUBRESFEATURES2D_MOVIE detects subresolution features in a series of images
%
%SYNOPSIS [movieInfo,emptyFrames,framesFailed,errFlag] = detectSubResFeatures2D_Movie(...
%     movieParam,detectionParam,saveResults)
%
%INPUT  movieParam    : Structure with fields
%           .imageDir     : Directory where images are stored
%           .candsDir     : Directory where cands (initial maxima) are stored.
%           .filenameBase : Filename base.
%           .firstImageNum: Numerical index of first image in movie.
%           .lastImageNum : Numerical index of last image in movie.
%           .digits4Enum  : Number of digits used to enumerate frames.
%       detectionParam: Structure with fields
%           .psfSigma     : Standard deviation of point spread function (in pixels).
%           .testAlpha    : Alpha-values for statistical tests. Optional.
%                           (See detectSubResFeatures2D for details).
%           .visual       : 1 if user wants to view results; 0 otherwise.
%                           Optional. Default: 0.
%           .doMMF        : 1 if user wants to do mixture-model fitting, 0
%                           otherwise. Optional. Default: 1.
%           .bitDepth     : Camera bit depth. Optional. Default: 14.
%       saveResults   : Structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                           Optional. Default: detectedFeatures.
%                       Whole structure optional.
%
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT movieInfo     : Array of length "movie length" of structures
%                       containing the fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       emptyFrames   : Array indicating frames where no features were
%                       detected.
%       framesFailed  : Array indicating frames where detection has failed.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, July 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movieInfo = [];
emptyFrames = [];
framesFailed = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--detectSubResFeatures2D_Movie: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get movie parameters
imageDir = movieParam.imageDir;
candsDir = movieParam.candsDir;
filenameBase = movieParam.filenameBase;
firstImageNum = movieParam.firstImageNum;
lastImageNum = movieParam.lastImageNum;
digits4Enum = movieParam.digits4Enum;

%get PSF sigma
psfSigma = detectionParam.psfSigma;

%get statistical test alpha values
if ~isfield(detectionParam,'testAlpha') || isempty(detectionParam.testAlpha)
    testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0.05);
else
    testAlpha = detectionParam.testAlpha;
end

%get visualization option
if ~isfield(detectionParam,'visual') || isempty(detectionParam.visual)
    visual = 0;
else
    visual = detectionParam.visual;
end

%check whether to do MMF
if ~isfield(detectionParam,'doMMF') || isempty(detectionParam.doMMF)
    doMMF = 1;
else
    doMMF = detectionParam.doMMF;
end

%get camera bit depth
if ~isfield(detectionParam,'bitDepth') || isempty(detectionParam.bitDepth)
    bitDepth = 14;
else
    bitDepth = detectionParam.bitDepth;
end

%determine where to save results
if nargin < 3 || isempty(saveResults) %if nothing was input
    saveResDir = pwd;
    saveResFile = 'detectedFeatures';
else
    if ~isfield(saveResults,'dir') || isempty(saveResults.dir)
        saveResDir = pwd;
    else
        saveResDir = saveResults.dir;
    end
    if ~isfield(saveResults,'filename') || isempty(saveResults.filename)
        saveResFile = 'detectedFeatures';
    else
        saveResFile = saveResults.filename;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign leading zeros in numerical index of images
switch digits4Enum
    case 4
        leadingZeros(1).value = '000';
        leadingZeros(2).value = '00';
        leadingZeros(3).value = '0';
        leadingZeros(4).value = '';
    case 3
        leadingZeros(1).value = '00';
        leadingZeros(2).value = '0';
        leadingZeros(3).value = '';
    case 2
        leadingZeros(1).value = '0';
        leadingZeros(2).value = '';
    case 1
        leadingZeros(1).value = '';
end

movieInfo(lastImageNum).xCoord = [];
movieInfo(lastImageNum).yCoord = [];
movieInfo(lastImageNum).amp = [];

%define variable to monitor progress
iProgress = 0;
progressText(0);

%go through the images

for i=min(9999,lastImageNum):-1:max(1000,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(4).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(4).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,...
            testAlpha,visual,doMMF,bitDepth);

        %save results
        movieInfo(i) = featuresInfo;

        %check whether frame is empty
        if isempty(featuresInfo.xCoord)
            emptyFrames = [emptyFrames; i];
        end

    catch %if detection fails

        %label frame as empty
        emptyFrames = [emptyFrames; i];

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];
        
    end

    %display progress on screen
    iProgress = iProgress + 1;
    progressText(iProgress/(lastImageNum-firstImageNum+1));

end

for i=min(999,lastImageNum):-1:max(100,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(3).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(3).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,...
            testAlpha,visual,doMMF,bitDepth);

        %save results
        movieInfo(i) = featuresInfo;

        %check whether frame is empty
        if isempty(featuresInfo.xCoord)
            emptyFrames = [emptyFrames; i];
        end

    catch %if detection fails

        %label frame as empty
        emptyFrames = [emptyFrames; i];

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

    %display progress on screen
    iProgress = iProgress + 1;
    progressText(iProgress/(lastImageNum-firstImageNum+1));

end

for i=min(99,lastImageNum):-1:max(10,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(2).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(2).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,...
            testAlpha,visual,doMMF,bitDepth);

        %save results
        movieInfo(i) = featuresInfo;

        %check whether frame is empty
        if isempty(featuresInfo.xCoord)
            emptyFrames = [emptyFrames; i];
        end

    catch %if detection fails

        %label frame as empty
        emptyFrames = [emptyFrames; i];

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

    %display progress on screen
    iProgress = iProgress + 1;
    progressText(iProgress/(lastImageNum-firstImageNum+1));

end

for i=min(9,lastImageNum):-1:max(1,firstImageNum)

    %get image
    image = imread([imageDir filenameBase leadingZeros(1).value num2str(i) '.tif']);

    %get cands
    eval(['load ' candsDir 'cands' leadingZeros(1).value num2str(i) ';'])

    try %try to detect features in this frame

        %fit with mixture-models
        featuresInfo = detectSubResFeatures2D(image,cands,psfSigma,...
            testAlpha,visual,doMMF,bitDepth);

        %save results
        movieInfo(i) = featuresInfo;

        %check whether frame is empty
        if isempty(featuresInfo.xCoord)
            emptyFrames = [emptyFrames; i];
        end

    catch %if detection fails

        %label frame as empty
        emptyFrames = [emptyFrames; i];

        %add this frame to the array of frames with failed detection
        framesFailed = [framesFailed; i];

    end

    %display progress on screen
    iProgress = iProgress + 1;
    progressText(iProgress/(lastImageNum-firstImageNum+1));

end

%save results
save([saveResDir filesep saveResFile],'movieParam','detectionParam',...
    'movieInfo','emptyFrames','framesFailed');


%%%%% ~~ the end ~~ %%%%%
