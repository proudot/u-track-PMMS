function [trackStats,statsRelChange,errFlag] = getTrackStats(trackedFeatureInfo,lenFrac,...
    timeWindow,problem2D,trackStatsOld)
%GETTRACKSTATS determines the statistical characeteristics of amplitude change and displacement in tracks over time
%
%SYNOPSIS [trackStats,errFlag] = getTrackStats(trackedFeatureInfo,lenFrac,...
%    timeWindow,problem2D,trackStatsOld)
%
%INPUT  trackedFeatureInfo:The positions and amplitudes of the tracked
%                          features. Number of rows = number of tracks, 
%                          while number of columns = 8*number of time 
%                          points. Each row consists of 
%                          [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                          in image coordinate system (coordinates in
%                          pixels). NaN is used to indicate time points 
%                          where the track does not exist.
%       lenFrac          : Minimum length of tracks used for statistical
%                          analysis, as a fraction of the total number of
%                          time points in movie.
%       timeWindow       : Time window of gap closing.
%                          Optional. Default: 1.
%       problem2D        : 1 if problem is 2D, 0 otherwise. Optional.
%                          Default: 0.
%       trackStatsOld    : trackStats (See output description) from previous 
%                          calculation. Optional. Default: [].
%
%OUTPUT trackStats       : Structure with fields:
%           .dispSqR        : timeWindow x 1 vector of r in the gamma
%                             distribution that describes the displacement
%                             of a feature between frames.
%           .dispSqTheta    : timeWindow x 1 vector of theta in the 
%                             gamma distribution that describes the 
%                             displacement of a feature between frames.
%           .ampDiffStd     : timeWindow x 1 vector of standard deviation of 
%                             the change in a feature's amplitude between
%                             frames.
%       statsRelChange   : Relative change in statistical parameters.
%       errFlag          : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackStats = [];
statsRelChange = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--getTrackStats: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%assign defaults
timeWindow_def = 1;

%check timeWindow
if nargin < 3 || isempty(timeWindow)
    timeWindow = timeWindow_def;
end

if nargin < 4 || isempty(problem2D)
    problem2D = 0;
end

if nargin < 5 || isempty(trackStatsOld)
    trackStatsOld = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of time points and tracks in movie
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/8;

%find tracks that contain more than (lenFrac*numTimePoints) time points
trackLength = zeros(numTracks,1);
for i=1:numTracks
    trackLength(i) = length(find(~isnan(trackedFeatureInfo(i,:))))/8;
end
goodTracks = find(trackLength>=lenFrac*numTimePoints);

%retain the lengths of only those tracks that are long enough
trackLength = trackLength(goodTracks);

%set curve fitting options
options = optimset('Display','off');

if ~isempty(goodTracks) %if there are tracks to use ...

    %get the x,y-coordinates and amplitudes over time of these good tracks
    xCoord = trackedFeatureInfo(goodTracks,1:8:end);
    yCoord = trackedFeatureInfo(goodTracks,2:8:end);
    zCoord = trackedFeatureInfo(goodTracks,3:8:end);
    amplitude = trackedFeatureInfo(goodTracks,4:8:end);

    %reserve memory for output vectors
    dispSqR = zeros(timeWindow,1);
    dispSqTheta = zeros(timeWindow,1);
    ampDiffStd = zeros(timeWindow,1);
    
    for i=1:timeWindow

        %get the squared displacement between time points
        dispSq = (xCoord(:,i+1:end) - xCoord(:,1:end-i)).^2 + ...
            (yCoord(:,i+1:end) - yCoord(:,1:end-i)).^2 + ...
            (zCoord(:,i+1:end) - zCoord(:,1:end-i)).^2;

        %if the problem is 2D, the squared displacements have an exponential
        %distribution, which is a special case of the gamma distribution. I
        %treat the exponential distribution explicitly in order to make the
        %tracking more reliable and robust.
        if problem2D

            %obtain the histogram of the squared displacement
            [n,x] = optimalHistogram(dispSq(:));

            %fit the histogram with an exponential function
            expParam = lsqcurvefit(@expFun,[100 1]',x,n,[],[],options);

            %assign gamma distribution parameters
            dispSqR(i) = 1;
            dispSqTheta(i) = expParam(2);

        else %if problem is not 2D

            %calculate the mean and standard deviation of the squared
            %displacement
            dispSqMean = sum(trackLength.*nanmean(dispSq,2))/sum(trackLength(:));
            dispSqVar = sum(trackLength.*nanvar(dispSq,[],2))./sum(trackLength(:));

            %calculate the gamma distribution parameters
            %avoid zeros because they cause problems with the logarithm and
            %gamma functions
            dispSqR(i) = max(dispSqMean^2/dispSqVar,1e-300); %gamma(realmin)=Inf, hence use 1e-300
            dispSqTheta(i) = max(dispSqMean/dispSqVar,realmin);

        end %(if problem2D)

        %get the amplitude change between time points
        ampDiff = amplitude(:,i+1:end) - amplitude(:,1:end-i);

        %get the standard deviation of amplitude change
        ampDiffStd(i) = max(nanstd(ampDiff(:)),realmin);

    end %(for i=1:timeWindow)

    %save output in structure
    trackStats.dispSqR = dispSqR;
    trackStats.dispSqTheta = dispSqTheta;
    trackStats.ampDiffStd = ampDiffStd;

    %get the maximum relative change in parameters, if the old parameters
    %are supplied
    if ~isempty(trackStatsOld)
        oldParam = [trackStatsOld.dispSqR; trackStatsOld.dispSqTheta; trackStatsOld.ampDiffStd];
        statsRelChange = max(abs(([dispSqR; dispSqTheta; ampDiffStd] - oldParam)...
            ./oldParam));
    end

else %if there aren't any ...

    errFlag = 1;
    return

end %(if ~isempty(goodTracks) ... else ...)


%%%%% ~~ the end ~~ %%%%%

