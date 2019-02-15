function [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal,movieInfo)
%CONVSTRUCT2MATNOMS converts tracks from structure format to matrix format, provided there are NO merges/splits.
%
%SYNPOSIS [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal,movieInfo)
%
%INPUT  tracksFinal: Output of trackCloseGapsKalman, when run with
%                    gapCloseParam.mergeSplit = 0.
%       trackedFeatureInfo, trackedFeatureIndx: Output of trackWithGapClosing.
%
%Khuloud Jaqaman, February 2008
%Kathryn Applegate, August 2009 - added movieInfo input to make larger 
%matrix from full movie, since if fewer frames than the total are tracked,
%numTimePoints will be smaller if we only look at the max frame number
%used

%% conversion

%get number of tracks
numTracks = length(tracksFinal);

%get number of time points
if nargin<2 || isempty(movieInfo)
    tmp = vertcat(tracksFinal.seqOfEvents);
    numTimePoints = max(tmp(:,1));
else
    numTimePoints=length(movieInfo);
end

%reserve memory for matrix of tracks
trackedFeatureInfo = NaN(numTracks,8*numTimePoints);

%put tracks in matrix
for iTrack = 1 : numTracks
    startTime = tracksFinal(iTrack).seqOfEvents(1,1);
    endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
    trackedFeatureInfo(iTrack,8*(startTime-1)+1:8*endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(1,:);
end

if nargout == 2
    
    %reserve memory for matrix of feature indices
    trackedFeatureIndx = zeros(numTracks,numTimePoints);
    
    %put indices in matrix
    for iTrack = 1 : numTracks
        startTime = tracksFinal(iTrack).seqOfEvents(1,1);
        endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
        trackedFeatureIndx(iTrack,startTime:endTime) = ...
            tracksFinal(iTrack).tracksFeatIndxCG(1,:);
    end
    
end

%% ~~~ the end ~~~