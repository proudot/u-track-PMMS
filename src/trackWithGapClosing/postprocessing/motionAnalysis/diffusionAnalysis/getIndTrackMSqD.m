function [mSqDisp,mSqDMatrix,mSqDErrMatrix,errFlag] = ...
    getIndTrackMSqD(trackedFeatureInfo,maxLag)
%GETINDTRACKMSQD calculates the mean squared displacement of the input tracks
%
%SYNOPSIS [mSqDisp,mSqDMatrix,mSqDErrMatrix,errFlag] = ...
%    getIndTrackMSqD(trackedFeatureInfo,maxLag)
%
%INPUT  trackedFeatureInfo: Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...].
%                           NaN is used to indicate time points where the
%                           track does not exist.
%       maxLag            : Maximum lag for mean square displacement
%                           calculation.
%                           Optiona. Default: track length/4.
%
%OUTPUT mSqDisp           : Structure array of length equal to number of
%                           input tracks with field "values" = output of
%                           the function meanSquaredDisplacement.
%       mSqDMatrix        : Matrix where the mean square displacement of
%                           each track is placed in a column. NaN is added
%                           to the end of each column to equalize column
%                           lengths.
%       mSqDErrMatrix     : Matrix where the error in the mean square
%                           displacement of each track is placed in a column.
%                           NaN also added for length equalization.
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mSqDisp = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getIndTrackMSqD: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

if nargin < 2
    maxLag = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mean squared displacement calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of tracks and time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/8;

%reserve memory
mSqDisp = repmat(struct('values',[]),numTracks,1);
timeLagMax = zeros(numTracks,1);

%go over all tracks
for j=numTracks:-1:1
    
    %get positions over time
    coordinates = [trackedFeatureInfo(j,1:8:end)' ...
        trackedFeatureInfo(j,2:8:end)' trackedFeatureInfo(j,3:8:end)'];

    %get corresponding variance-covariance matrices
    for i=1:numTimePoints
        standardDevs = [trackedFeatureInfo(j,5:8:end)' ...
            trackedFeatureInfo(j,6:8:end)' trackedFeatureInfo(j,7:8:end)'];
    end

    %get the mean squared displacement
    trackMoments = calcTrackMoments(coordinates,standardDevs,2,maxLag);
    mSqDisp(j).values = trackMoments.momentValues;

    %get the maximum time lag
    timeLagMax(j) = size(trackMoments.momentValues,1);

end

%collect all mean square displacements and their errors in two big matrices

%get the maximum time lag among all tracks
timeLagMaxMax = max(timeLagMax);

%initialize the two big matrices
mSqDMatrix = NaN*ones(timeLagMaxMax,numTracks);
mSqDErrMatrix = NaN*ones(timeLagMaxMax,numTracks);

%fill the matrices
for j=1:numTracks
    mSqDMatrix(1:timeLagMax(j),j) = mSqDisp(j).values(:,1); %mean square displacement
    mSqDErrMatrix(1:timeLagMax(j),j) = mSqDisp(j).values(:,2); %its standard deviation
end


%%%%% ~~ the end ~~ %%%%%
