function trackSEL = getTrackSEL(trackedFeatureInfo)
%GETTRACKSEL outputs track start times, end times and lifetimes
%
%SYNOPSIS trackSEL = getTrackSEL(trackedFeatureInfo);
%
%INPUT  trackedFeatureInfo: Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%
%OUTPUT trackSEL          : An array with 3 columns and number of rows equal
%                           to number of tracks. 1st column indicates track 
%                           start times, 2nd column indicates track end
%                           times and 3rd column indicates track lifetimes.
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackSEL = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getTrackSEL: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track information extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make new matrix which contains only one column per time point
trackedFeatureInfo = trackedFeatureInfo(:,1:8:end);

%get number of tracks
numTracks = size(trackedFeatureInfo,1);

%alocate memory for output
trackSEL = zeros(numTracks,3);

%find track start times
for i=1:numTracks
    trackSEL(i,1) = find(~isnan(trackedFeatureInfo(i,:)),1,'first');
end

%find track end times
for i=1:numTracks
    trackSEL(i,2) = find(~isnan(trackedFeatureInfo(i,:)),1,'last');
end

%calculate track lifetimes
trackSEL(:,3) = trackSEL(:,2) - trackSEL(:,1) + 1;


%%%%% ~~ the end ~~ %%%%%

