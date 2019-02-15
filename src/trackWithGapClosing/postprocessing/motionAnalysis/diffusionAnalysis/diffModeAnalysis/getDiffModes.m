function [modeParam,numMode,modeParamControl,numModeControl] = ...
    getDiffModes(tracksFinal,minLength,alpha,showPlot,maxNumMode,...
    binStrategy,plotName,subSampSize,doControl)
%GETDIFFMODES determines number of diffusion modes and their parameters from distribution of frame-to-frame displacements
%
%SYNOPSIS [modeParam,numMode,modeParamControl,numModeControl] = ...
%    getDiffModes(tracksFinal,minLength,alpha,showPlot,maxNumMode,...
%    binStrategy,plotName,subSampSize)
%
%INPUT  tracksFinal : Output of trackCloseGapsKalman.
%                     Optional. If not input, GUI will be called to get
%                     tracks. NOTE: This only works if tracks variable is
%                     called tracksFinal inside file.
%       minLength   : Minimum length of a track to be included in analysis.
%                     Optional. Default: 5 frames.
%       alpha       : Alpha-value for the statistical test to determine
%                     number of modes.
%                     Optional. Default: 0.01.
%       showPlot    : 0 to not plot anything.
%                     1 to plot the histogram and fitted exponentials.
%                     Optional. Default: 1.
%       maxNumMode  : Upper limit on the number of modes.
%                     Optional. Default: 10.            
%       binStrategy : Binning strategy for calculating the cumulative
%                     histogram. 1 for using "histogram" and 2 for using
%                     the data directly.
%                     Optional. Default: 2.
%       plotName    : The title of the plotted figure.
%                     Optional. Default: 'Figure'.
%       subSampSize : Size of subsample to use in mode decomposition. In
%                     this case, the original data are subsampled many
%                     times, each time with the specified subSampSize. The
%                     output is then the diffusion mode decomposition
%                     result for each of the subsamples. Enter [] if
%                     no sub-sampling.
%                     Optional. Default: [].
%       doControl   : 1 in order to do mono-exponential control, 0
%                     otherwise. 
%                     Optional. Default: 1.
%                     
%OUTPUT modeParam   : Matrix with number of rows equal to number of modes
%                     and 4 columns:
%                     Column 1: diffusion coefficient of each mode.
%                     Column 2: fraction of contribution of each mode.
%                     Column 3: std of each diffusion coefficient.
%                     Column 4: std of each fraction of contribution.
%                     In case of subsampling, the 3rd dimension refers to
%                     the results of each subsample.
%       numMode     : Number of diffusion modes. In case of subsampling,
%                     this is an array with a value for each subsample.
%       modeParamControl: The same diffusion mode decomposition but
%                     performed on synthetic data representing a truly
%                     mono-exponential distribution. The synthetic data
%                     have the same size and mean as the input data.
%       numModeControl: Number of diffusion modes in control data.
%
%REMARKS Code is written for 2D case only, but can be generalized to 3D.
%
%Khuloud Jaqaman, June 2012

%% Input

%get tracks if not input
%this only works if tracks are called tracksFinal inside the .mat file
if nargin < 1 || isempty(tracksFinal)
    [fName,dirName] = uigetfile('*.mat','Please choose .mat file of tracks.');
    load(fullfile(dirName,fName));
end

%check rest of input

if nargin < 2 || isempty(minLength)
    minLength = 5;
end

if nargin < 3 || isempty(alpha)
    alpha = 0.01;
end

if nargin < 4 || isempty(showPlot)
    showPlot = 1;
end

if nargin < 5 || isempty(maxNumMode)
    maxNumMode = 10;
end

if nargin < 6 || isempty(binStrategy)
    binStrategy = 2;
end

if nargin < 7 || isempty(plotName)
    plotName = 'Figure';
end

if nargin < 8 || isempty(subSampSize)
    subSampSize = [];
end
if isempty(subSampSize)
    numSamp = 1;
else
    numSamp = 10;
    showPlot = 0;
end

if nargin < 9 || isempty(doControl)
    doControl = 1;
end

%% Square displacement distribution

%keep only tracks of the right length
criteria.lifeTime.min = minLength;
indxKeep = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indxKeep);

%get number of tracks
numTracks = length(tracksFinal);

%get distribution of square displacements
%divide tracks into chunks of 5000 so as not to run out of memory
%also get localization precision
indxFirst = 1;
f2fdispSqAll = [];
xCoordStdAll = [];
yCoordStdAll = [];
while indxFirst < numTracks

    %get chunk of tracks for this iteration
    indxLast = min(indxFirst+5000-1,numTracks);
    tracksTmp = tracksFinal(indxFirst:indxLast);
    tracksMat = convStruct2MatIgnoreMS(tracksTmp,1);

    xCoord = tracksMat(:,1:8:end);
    yCoord = tracksMat(:,2:8:end);
    f2fdispSqTmp = diff(xCoord,1,2).^2 + diff(yCoord,1,2).^2;
    f2fdispSqTmp = f2fdispSqTmp(~isnan(f2fdispSqTmp));
    f2fdispSqAll = [f2fdispSqAll; f2fdispSqTmp]; %#ok<AGROW>

    xCoordStd = tracksMat(:,5:8:end);
    yCoordStd = tracksMat(:,6:8:end);
    xCoordStdAll = [xCoordStdAll; xCoordStd(~isnan(xCoordStd))]; %#ok<AGROW>
    yCoordStdAll = [yCoordStdAll; yCoordStd(~isnan(yCoordStd))]; %#ok<AGROW>
    
    indxFirst = indxLast + 1;

end
meanPosVar = mean([xCoordStdAll;yCoordStdAll].^2);
dataSize = length(f2fdispSqAll);

%generate a mono-exponenetial distribution with the same mean as a control
f2fdispSqControl = exprnd(mean(f2fdispSqAll),dataSize,1);

%patch-up for for-loop later
if numSamp == 1
    subSampSize = dataSize;
end

%% Exponential fits

%initialize memory
numMode = zeros(numSamp,1);
numModeControl = zeros(numSamp,1);
expParam = zeros(maxNumMode,6,numSamp);
expParamControl = zeros(maxNumMode,6,numSamp);

%go over all subsamples
for iSamp = 1 : numSamp
    
    %get subsample
    indxSamp = randsample(dataSize,subSampSize);
    
    %fit data
    expParamTmp = fitHistWithExponentialsN(f2fdispSqAll(indxSamp),...
        alpha,showPlot,maxNumMode,binStrategy,[],plotName,4*meanPosVar);
    numMode(iSamp) = size(expParamTmp,1);
    expParam(1:numMode(iSamp),:,iSamp) = expParamTmp;
    
    %fit control
    if doControl
        expParamTmp = fitHistWithExponentialsN(f2fdispSqControl(indxSamp),...
            alpha,showPlot,maxNumMode,binStrategy,[],plotName,4*meanPosVar);
        numModeControl(iSamp) = size(expParamTmp,1);
        expParamControl(1:numModeControl(iSamp),:,iSamp) = expParamTmp;
    end
    
end

%remove unnecessary rows
numModeMaxData = max(numMode);
expParam = expParam(1:numModeMaxData,:,:);
expParam(expParam==0) = NaN;
if doControl
    numModeMaxControl = max(numModeControl);
    expParamControl = expParamControl(1:numModeMaxControl,:,:);
    expParamControl(expParamControl==0) = NaN;
end

%% Diffusion modes

%data

%calculate diffusion coefficient of each mode
%formula: mu (i.e. mean of exponential) = 4*diffCoef + 4*posStd^2
diffCoef = expParam(:,1,:)/4 - meanPosVar;

%calculate the diffusion coefficient std -- under-estimate
diffCoefStd = expParam(:,3,:)/4;

%calculate fraction of contribution of each mode
fracContr = expParam(:,2,:)./repmat(nansum(expParam(:,2,:),1),[numModeMaxData 1 1]);

%calculate corresponding std -- under-estimate
fracContrStd = expParam(:,4,:)./repmat(nansum(expParam(:,2,:),1),[numModeMaxData 1 1]);

%output
modeParam = cat(2,diffCoef,fracContr,diffCoefStd,fracContrStd);

%control
if doControl
    diffCoef = expParamControl(:,1,:)/4 - meanPosVar;
    diffCoefStd = expParamControl(:,3,:)/4;
    fracContr = expParamControl(:,2,:)./repmat(nansum(expParamControl(:,2,:),1),[numModeMaxControl 1 1]);
    fracContrStd = expParamControl(:,4,:)./repmat(nansum(expParamControl(:,2,:),1),[numModeMaxControl 1 1]);
    modeParamControl = cat(2,diffCoef,fracContr,diffCoefStd,fracContrStd);
else
    modeParamControl = [];
    numModeControl = [];
end

%% ~~~ the end ~~~



