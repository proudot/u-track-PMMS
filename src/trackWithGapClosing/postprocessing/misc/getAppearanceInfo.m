function [appearanceInfo,errFlag] = getAppearanceInfo(tracksFinal,probDim)
%GETAPPEARANCEINFO extracts the appearance position of particles in a time-lapse sequence
%
%SYNOPSIS [appearanceInfo,errFlag] = getAppearanceInfo(tracksFinal,probDim)
%
%INPUT  trackFinal: Output of trackCloseGapsKalman.
%       probDim   : Problem dimensionality (2 or 3).
%                   Optional. Default: 2.
%
%OUTPUT appearanceInfo: A variable with the same fields as movieInfo
%                       (output of detectSubResFeatures2D_StandAlone) but
%                       with only the information on particle appearances.
%
%Khuloud Jaqaman, March 2011

%% Output

appearanceInfo = [];
errFlag = 0;

%% Input

%check number of input arguments
if nargin < 1
    disp('--getAppearanceInfo: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%check problem dimensionality
if nargin < 2 || isempty(probDim)
    probDim = 2;
else
    if probDim ~= 2 && probDim ~= 3
        disp('--getAppearanceInfo: Problem dimensionality can be only 2 or 3.');
        errFlag = 1;
        return
    end
end

%% Extraction of appearance information

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%find number of frames in movie
numFrames = max(trackSEL(:,2));

%initialize output variable
if probDim == 2
    appearanceInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numFrames,1);
else
    appearanceInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),numFrames,1);
end

%go over all frames
for iFrame = 1 : numFrames
    
    %find all tracks that start in this frame
    trackIndx = find(trackSEL(:,1)==iFrame);
    numTracks = length(trackIndx);
    
    %get their initial position and amplitude
    [xCoord,yCoord,zCoord,amp] = deal(zeros(numTracks,2));
    for iTrack = 1 : numTracks
        trackInfo = tracksFinal(trackIndx(iTrack)).tracksCoordAmpCG(1,1:8);
        xCoord(iTrack,:) = [trackInfo(1) trackInfo(5)];
        yCoord(iTrack,:) = [trackInfo(2) trackInfo(6)];
        zCoord(iTrack,:) = [trackInfo(3) trackInfo(7)];
        amp(iTrack,:)    = [trackInfo(4) trackInfo(8)];
    end

    %store information in output variable
    appearanceInfo(iFrame).xCoord = xCoord;
    appearanceInfo(iFrame).yCoord = yCoord;
    if probDim == 3
        appearanceInfo(iFrame).zCoord = zCoord;
    end
    appearanceInfo(iFrame).amp = amp;
    
end

%% ~~~ the end ~~~
