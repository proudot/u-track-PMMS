function [meanSqDisp,errFlag] = getAllTracksMSqD(tracks,maxLag,ensemble)
%GETALLTRACKSMSQD calculates the mean squared displacement of the input tracks all together
%
%SYNOPSIS [meanSqDisp,errFlag] = getAllTracksMSqD(tracks,maxLag,ensemble)
%
%INPUT  tracks            : Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%       maxlag            : Maximum lag for calculating mean squared
%                           displacement.
%       ensemble          : 1 to use ensemble average, 0 to use time
%                           average (default).
%
%OUTPUT meanSqDisp        : maxLag - by - 3 array. Row i is for lag i.
%                           1st column: mean of square displacements.
%                           2nd column: std of square displacements.
%                           3rd column: number of observations used in
%                               calculation.
%
%REMARKS Ensemble averaging is standard ensemble averaging. Time averaging
%        is taken from Saxton, Biophys. J., 1997, using overlapping time
%        intervals.
%
%Khuloud Jaqaman, November 2007

%% Input

%check whether correct number of input arguments was used
errFlag = 0;
if nargin < 2
    disp('--getAllTracksMSqD: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

if nargin < 3 || isempty(ensemble)
    ensemble = 0;
end

%% Mean square displacement calculation

%put all x-coordinates in one matrix, all y-coordinates in one matrix, etc.
xCoord = tracks(:,1:8:end);
yCoord = tracks(:,2:8:end);
zCoord = tracks(:,3:8:end);
clear tracks

%allocate memory for output
meanSqDisp = NaN(maxLag,3);

%go over all lags ...
for iLag = 1 : maxLag
    
    %construct the displacement vectors
    xCoordDelta2 = (xCoord(:,iLag+1:end) - xCoord(:,1:end-iLag)).^2;
    yCoordDelta2 = (yCoord(:,iLag+1:end) - yCoord(:,1:end-iLag)).^2;
    zCoordDelta2 = (zCoord(:,iLag+1:end) - zCoord(:,1:end-iLag)).^2;
    
    %calculate the square displacements
    squaredDisp = xCoordDelta2 + yCoordDelta2 + zCoordDelta2;
    if ensemble
        squaredDisp = squaredDisp(:,1);
    else
        squaredDisp = squaredDisp(:);
    end
    squaredDisp = squaredDisp(~isnan(squaredDisp));
    
    %calculate the mean square displacement, its standard deviation, and
    %the number of observations used in the calculation
    meanSqDisp(iLag,:) = [mean(squaredDisp) std(squaredDisp) ...
        length(squaredDisp)];
    
end

%% %%%%% ~~ the end ~~ %%%%%
