function [movieInfoGT,tracksGT,errFlag] = simEB1Images(imSize,pixelSize,...
    mtDensity,SNR,ampAboveBG,samplingRate,totalTime,mtSimParam,ebCometParam,saveInfo)
%simEB1Images generates images of EB1 comets on MTs undergoing dynamic instability
%
%[movieInfoGT,tracksGT,errFlag] = simEB1Images(imSize,pixelSize,...
%    mtDensity,SNR,ampAboveBG,samplingRate,totalTime,mtSimParam,ebCometParam,saveInfo)
%
%INPUT  imSize      : Image size in x = y, in pixels.
%       pixelSize   : Pixel size, in micrometers.
%       mtDensity   : Density of MTs, in number of MTs / pixel.
%       SNR         : Image signal-to-noise ratio.
%       ampAboveBG  : EB1 comet amplitude above background.
%       samplingRate: Image sampling rate, in Hz.
%       totalTime   : Total time of simulation [s].
%       mtSimPara   : Parameters for simulating MT dynamics. Same as
%                     modelParam in mtGammaTdSd.
%       ebCometParam: Sigmas of Gaussian approximating EB comet, along the
%                     long and short axes.
%       saveInfo    : Information for saving the EB1 images. Structure with
%                     fields:
%           .filenameBase: Names of image tif files.
%           .dir2save    : Name of directory where results are to be saved.

%
%OUTPUT movieInfoGT : Ground truth positions in movieInfo format.
%       tracksGT    : Ground truth tracks in tracksFinal format.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%       The code also generates a series of images that get named and
%       stored as indicated in saveInfo.
%
%Khuloud Jaqaman, June 2011

%% Output

movieInfoGT = [];
tracksGT = [];
errFlag = 0;

%% Input

if nargin < nargin('simEB1Images')
    disp('--simEB1Images: Wrong number of input variables!');
    errFlag = 1;
    return
end

%% Simulation

%calculate time between frames and total number of frames
frameInterval = 1 / samplingRate;
numFrames = ceil(totalTime / frameInterval);
totalTimeMod = (numFrames+100) * frameInterval;

%get number of MTs to simulate
numMTs = round(mtDensity * imSize^2);

%assign the MTs random initial lengths, between 0 and imSize
mtLength = rand(numMTs,1) * imSize * pixelSize;

%also assign MTs random orientations, between 0 and pi
orientEB1 = rand(numMTs,1) * pi/2;

%go over all MTs and get their lengths and EB1 positions
positionMT = NaN(numFrames,2,numMTs);
positionEB1 = NaN(numFrames,2,numMTs);
for iMT = 1 : numMTs
    
    %simulate the MT dynamics
    [simTraj] = mtGammaTdSd(mtSimParam,mtLength(iMT),totalTimeMod,[0 min(mtLength(iMT)+100*pixelSize,imSize*pixelSize)]);
    
    %sample the MT length at the "imaging" frames
    [trajSamp] = sampleTraj(simTraj,frameInterval);
    trajSamp = trajSamp(:,2);
    
    %assign +1 to growth phases and -1 to shrinkage phases
    phaseMT = [1; sign(diff(trajSamp))];
    
    %from the oriention, calculate the EB1 positions
    xCoordTime = trajSamp * cos(orientEB1(iMT));
    yCoordTime = trajSamp * sin(orientEB1(iMT));
    
    %put this MT's values in the big matrix, after converting to pixels
    positionMT(:,1,iMT) = xCoordTime(101:numFrames+100) / pixelSize;
    positionMT(:,2,iMT) = yCoordTime(101:numFrames+100) / pixelSize;
    
    %now generate the EB1 positions, which don't exist for shrinking MTs
    positionEB1(:,:,iMT) = positionMT(:,:,iMT);
    positionEB1(phaseMT(101:numFrames+100)==-1,:,iMT) = NaN;
    
end

%% Image generation

%reserve memory for movieInfoGT
movieInfoGT = repmat(struct('xCoord',[],'yCoord',[],'orient',[],'amp',[]),numFrames,1);

%calculate number of digits for image enumeration
numDigits = num2str(ceil(log10(numFrames))+1);

%estimate noise standard deviation from SNR
noiseStd = ampAboveBG / SNR;

%go over all time points
for iFrame = 1 : numFrames
    
    %generate movieInfoGT, which contains for now positions and orientations
    %(no amplitude information)
    xCoordFrame = squeeze(positionEB1(iFrame,1,:));
    yCoordFrame = squeeze(positionEB1(iFrame,2,:));
    orientFrame = orientEB1(~isnan(xCoordFrame));
    xCoordFrame = xCoordFrame(~isnan(xCoordFrame));
    yCoordFrame = yCoordFrame(~isnan(yCoordFrame));
    movieInfoGT(iFrame).xCoord = [xCoordFrame 0.5*ones(size(xCoordFrame))];
    movieInfoGT(iFrame).yCoord = [yCoordFrame 0.5*ones(size(xCoordFrame))];
    movieInfoGT(iFrame).orient = [orientFrame zeros(size(xCoordFrame))];
    movieInfoGT(iFrame).amp = [ones(size(xCoordFrame)) 0.1*ones(size(xCoordFrame))];
    
    %get number of features in this frame
    numFeat = length(xCoordFrame);
    
    %go over these features and generate the image, initially with no noise
    imageEB1 = zeros(imSize,imSize);
    for iFeat = 1 : numFeat
        [xRange,yRange,nzIdx] = anisoGaussian2DSupport(xCoordFrame(iFeat),yCoordFrame(iFeat),ebCometParam(1),ebCometParam(2),orientFrame(iFeat),4,[imSize imSize]);
        if ~isempty(nzIdx)
            F = anisoGaussian2D(xCoordFrame(iFeat),yCoordFrame(iFeat),ampAboveBG,ebCometParam(1),ebCometParam(2),orientFrame(iFeat),xRange,yRange,nzIdx);
            myImage = zeros(numel(yRange),numel(xRange));
            myImage(nzIdx) = F;
            imageEB1(yRange,xRange) = imageEB1(yRange,xRange) + myImage;
        end
    end
    
    %now add noise and make sure that all intensities are above zero
    noiseImage = randn(imSize,imSize) * noiseStd;
    minNoiseVal = min(noiseImage(:));
    noiseImage = noiseImage - minNoiseVal + 1;
    imageEB1 = imageEB1 + noiseImage;
    
    %save image where specified
    imageFileName = [saveInfo.filenameBase '_' num2str(iFrame,['%0' numDigits 'i']) '.tif'];
    imwrite(uint16(round(imageEB1)),fullfile(saveInfo.dir2save,imageFileName),'tif','Compression','None');
    
end

%% EB1 tracks in the form of tracksFinal

%reserve memory for tracksGT
%tracksFeatIndxCG is not filled for now
tracksGT = repmat(struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],'seqOfEvents',[]),numMTs,1);

%store the track information for each EB1 comet
for iMT = 1 : numMTs

    %get this track's position information
    tracksCoordAmpCG = NaN(1,numFrames*8);
    goodIndx = find(~isnan(positionEB1(:,1,iMT)));
    goodTimes = (find(~isnan(positionEB1(:,1,iMT)))-1) * 8;
    tracksCoordAmpCG(goodTimes+1) = positionEB1(goodIndx,1,iMT);
    tracksCoordAmpCG(goodTimes+2) = positionEB1(goodIndx,2,iMT);
    tracksCoordAmpCG(goodTimes+3) = 0;
    tracksCoordAmpCG(goodTimes+4) = 1;
    tracksCoordAmpCG(goodTimes+5) = 0.5;
    tracksCoordAmpCG(goodTimes+6) = 0.5;
    tracksCoordAmpCG(goodTimes+7) = 0;
    tracksCoordAmpCG(goodTimes+8) = 0.1;

    %find the start and end time of this track
    trackSEL = getTrackSEL(tracksCoordAmpCG);
    startColumn = (trackSEL(1)-1)*8 + 1;
    endColumn = trackSEL(2)*8;
    
    %modify tracksCoordAmpCG if necessary
    tracksCoordAmpCG = tracksCoordAmpCG(startColumn:endColumn);
    
    %store information for output
    tracksGT(iMT).tracksCoordAmpCG = tracksCoordAmpCG;
    tracksGT(iMT).seqOfEvents = [trackSEL(1) 1 1 NaN; trackSEL(2) 2 1 NaN];
    
end

