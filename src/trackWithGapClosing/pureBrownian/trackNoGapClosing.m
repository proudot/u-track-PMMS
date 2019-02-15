function [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
    trackNoGapClosing(movieInfo,linkCriteria,miscParam)
%TRACKNOGAPCLOSING links features between frames in a movie without closing gap
%
%SYNOPSIS function [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
%    trackNoGapClosing(movieInfo,linkCriteria,miscParam)
%
%INPUT  movieInfo       : Array of size equal to the number of time points
%                         in a movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinate of detected
%                            features [x dx] (in pixels).
%             .yCoord      : Image coorsinate system y-coordinate of detected
%                            features [y dy] (in pixels).
%             .amp         : Amplitudes of PSFs fitting detected features [a da].
%       linkCriteria    : Structure of linking parameters. Contains the
%                         fields:
%             .searchRadius: Maximum distance between two features in two
%                            consecutive time points that allows linking 
%                            them (in pixels).
%             .maxAmpRatio : Maximum ratio between the amplitudes of two
%                            features in two censecutive time points that 
%                            allows linking them.
%             .cutoffCProb : Cumulative probability of squared displacement
%                            or amplitude change beyond which links are not
%                            allowed.
%       miscParam       : Structure with the fields:
%             .tolerance   : Tolerance for changes in track statistics to
%                            stop iterating.
%             .lenFrac     : Minimum length of tracks used for statistical
%                            analysis, as a fraction of the total number of
%                            time points in movie.
%
%OUTPUT trackedFeatureNum: Connectivity matrix of features between time points.
%                          Rows indicate continuous tracks, while columns 
%                          indicate time points. A track that ends before the
%                          last time point is followed by zeros, and a track
%                          that starts at a time after the first time point
%                          is preceded by zeros. 
%       trackedFeatureInfo:The positions and amplitudes of the tracked
%                          features. Number of rows = number of tracks, 
%                          while number of columns = 6*number of time 
%                          points. Each row consists of 
%                          [x1 y1 a1 dx1 dy1 da1 x2 y2 a2 dx2 dy2 da2 ...]
%                          in image coordinate system (coordinates in
%                          pixels). NaN is used to indicate time points 
%                          where the track does not exist.
%       errFlag          : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm is currently for the special case of 2D, but in
%        principle it can be generalized to ND quite easily.
%
%MUST FIX!
%
%Khuloud Jaqaman, March 2006

%This function is a little outdated. Must be fixed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackedFeatureNum = [];
trackedFeatureInfo = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('trackNoGapClosing')
    disp('--trackNoGapClosing: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

searchRadius = linkCriteria.searchRadius;
maxAmpRatio = linkCriteria.maxAmpRatio;
cutoffCProb = linkCriteria.cutoffCProb;
tol = miscParam.tolerance;
lenFrac = miscParam.lenFrac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%link features between time points using a simple cost matrix
costMatParams = struct('searchRadius',searchRadius,'maxAmpRatio',maxAmpRatio);
[trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
    linkFeaturesTp2Tp(movieInfo,'costMatSimple',costMatParams);

%get track statistics
[dispSqLambdaT,ampDiffStdT,errFlag] = getTrackStats(trackedFeatureInfo,lenFrac);

%exit at this point if statistical analysis could not be performed
if errFlag
    disp('--trackNoGapClosing: getTrackStats failed. Only first iteration using simple cost matrix performed!');
    iterate = 0;
else %otherwise iterate
    iterate = 1;
    relParamChange = [1 1];
end

while iterate

    %update parameters
    dispSqLambda = dispSqLambdaT;
    ampDiffStd = ampDiffStdT;
    
    %get the maximum squared displacement that allows linking 2 features
    maxDispSq = expinv(cutoffCProb,dispSqLambda);

    %find the maximum squared amplitude change that allows linking 2 features
    maxAmpDiffSq = (norminv(cutoffCProb,0,ampDiffStd))^2;

    %link features between time points again using a cost matrix that uses the
    %distributions of amplitude change and squared displacement
    costMatParams = struct('dispSqLambda',dispSqLambda,'ampDiffStd',...
        ampDiffStd,'maxDispSq',maxDispSq,'maxAmpDiffSq',maxAmpDiffSq);
    [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
        linkFeaturesTp2Tp(movieInfo,'costMatLogL',costMatParams);

    %get track statistics
    [dispSqLambdaT,ampDiffStdT,errFlag] = getTrackStats(trackedFeatureInfo,lenFrac);

    %exit at this point if statistical analysis failed
    if errFlag
        disp('--trackNoGapClosing: getTrackStats failed. Stopping prematurely!');
        disp(['         Rel. param. change in last iteration: ' num2str(relParamChange)]);
        iterate = 0;
    else %otherwise check change in parameters and exit if smaller than tol
        relParamChange = [abs((dispSqLambdaT-dispSqLambda)/dispSqLambda) ...
            abs((ampDiffStdT-ampDiffStd)/ampDiffStd)];
        if min(relParamChange < tol)
            iterate = 0;
        end
    end

end %(while iterate)


%%%%% ~~ the end ~~ %%%%%
