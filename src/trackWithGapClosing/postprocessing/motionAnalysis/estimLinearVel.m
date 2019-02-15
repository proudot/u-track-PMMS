function speedDistr = estimLinearVel(tracks,kalmanInfo,minTrackLen,probDim)
%ESTIMLINEARVEL estimates the speeds in linear tracks
%
%SYNOPSIS velDistr = estimLinearVel(tracks,minTrackLen,probDim)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       kalmanInfo : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%
%OUTPUT speedDistr : Distribution of linear speed.
%
%Khuloud Jaqaman, December 2007

%% input

if nargin < 2 || isempty(tracks)
    disp('estimLinearVel: Missing input argument!');
    return
end

if nargin < 3 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 4 || isempty(probDim)
    probDim = 2;
end
if probDim == 3
    disp('estimLinearVel: Cannot deal with 3D yet');
    return
end 

%% preamble

%keep only linear tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
criteria.trackType = 1;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%% speed estimation

%initialization
speedDistr = [];

%go over all linear tracks ...
for iTrack = 1 : numTracks
    
    %construct matrix of linked features
    linkedFeatMat = tracks(iTrack).tracksFeatIndxCG;
    numSegments = size(linkedFeatMat,1);
    seqOfEvents = tracks(iTrack).seqOfEvents;
    linkedFeatMat = [zeros(numSegments,seqOfEvents(1,1)-1) ...
        linkedFeatMat zeros(numSegments,numFrames-seqOfEvents(end,1))];
    
    %get track's Kalman filter information
    trackKalmanInfo = getKalmanInfoLinearMotion2D(linkedFeatMat,kalmanInfo);
    
    %extract track's x and y velocities
    trackXVel = trackKalmanInfo(:,1:8:end);
    trackYVel = trackKalmanInfo(:,3:8:end);
    
    %calculate speed
    trackSpeed = sqrt( trackXVel.^2 + trackYVel.^2 );
    
    %retain only speeds of segments that are at least minTrackLen long
    segmentLen = zeros(numSegments,1);
    for iSegment = 1 : numSegments
        segmentLen(iSegment) = length(find(~isnan(trackSpeed(iSegment,:))));
    end
    trackSpeed = trackSpeed(segmentLen>=minTrackLen,:);
    
    %remove all NaNs
    trackSpeed = trackSpeed(~isnan(trackSpeed));
    
    %add speeds to overall distribution
    speedDistr = [speedDistr; trackSpeed(:)];
    
end

%% ~~~ the end ~~~
