function [trackedFeatKalmanInfo,errFlag] = getKalmanInfoLinearMotion2D(...
    trackedFeatureIndx,kalmanFilterInfo,indxList)
%GETKALMANINFOLINEARMOTION2D outputs the Kalman filter information for tracked features
%
%SYNOPSIS [trackedFeatKalmanInfo,errFlag] = getKalmanInfoLinearMotion2D(...
%    trackedFeatureIndx,kalmanFilterInfo,indxList)
%
%INPUT  trackedFeatureIndx: Connectivity matrix of features between time
%                           points. Rows indicate tracks, while columns
%                           indicate frames. A track that ends before the
%                           last frame is followed by zeros, and a track
%                           that starts at a time after the first frame
%                           is preceded by zeros. 
%       kalmanFilterInfo  : Kalman filter information as calculated in
%                           linkFeaturesKalman for the linear motion model.
%       indxList          : Row vector of indices of tracks whose Kalman
%                           filter information is to be extracted.
%                           Optional. Default: All tracks.
%
%OUTPUT trackedFeatKalmanInfo: The velocities, random element variance and
%                              propagation scheme along the tracks. Rows
%                              indicate tracks. # columns = 8 * # time
%                              points. Each row consists of
%                              [vx1 dvx1 vy1 dvy1 epsX1 epsY1 varEps1 scheme1 vx2 dvx2 vy2 dvy2 epsX2 epsY2 varEps2 scheme2 ...]
%
%       errFlag              : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackedFeatKalmanInfo = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--getKalmanInfoLinearMotion2D: Incorrect number of input arguments!');
    return
end

%get number of tracks and frames
[numTracks,numFrames] = size(trackedFeatureIndx);

%check indxList
if nargin < 3 || isempty(indxList)
    indxList = (1:numTracks);
else
    if size(indxList,1) ~= 1
        indxList = indxList';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalman information extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of specified tracks
numChosenTracks = length(indxList);

%reserve memory for output
trackedFeatKalmanInfo = NaN*ones(numChosenTracks,8*numFrames);

%go over all specified tracks
for iChosen = 1 : numChosenTracks

    %get track index
    iTrack = indxList(iChosen);

    %get connectivity vector of features making this track
    currTrack = trackedFeatureIndx(iTrack,:);

    %go over all frames
    for iFrame = 1 : numFrames

        %find index of feature making up the track in this frame (0
        %indicates that track does not exist in current frame)
        iFeature = currTrack(iFrame);

        %if the track exists in this frame
        if iFeature ~= 0

            %extract its Kalman filter information
            trackedFeatKalmanInfo(iChosen,8*(iFrame-1)+1:8*iFrame) = [...
                kalmanFilterInfo(iFrame).stateVec(iFeature,3) ... %x-component of velocity
                sqrt(kalmanFilterInfo(iFrame).stateCov(3,3,iFeature)) ... %its std
                kalmanFilterInfo(iFrame).stateVec(iFeature,4) ... %y-component of velocity
                sqrt(kalmanFilterInfo(iFrame).stateCov(4,4,iFeature)) ... %its std
                kalmanFilterInfo(iFrame).stateNoise(iFeature,1) ... %noise in x
                kalmanFilterInfo(iFrame).stateNoise(iFeature,2) ... %noise in y
                kalmanFilterInfo(iFrame).noiseVar(1,1,iFeature) ... %noise variance
                kalmanFilterInfo(iFrame).scheme(iFeature,2)]; %propagation scheme

        end %(if iFeature ~= 0)

    end %(for iFrame = 1 : numFrames)

end %(for iChosen = 1 : numChosenTracks)


%%%%% ~~ the end ~~ %%%%%
