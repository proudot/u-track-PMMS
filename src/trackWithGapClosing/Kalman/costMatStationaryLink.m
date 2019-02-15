function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
    errFlag] = costMatStationaryLink(movieInfo,kalmanFilterInfoFrame1,...
    costMatParam,nnDistFeatures,probDim,prevCost,featLifetime,...
    trackedFeatureIndx,currentFrame)
%costMatStationaryLink provides a cost matrix for linking stationary features
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
%     errFlag] = costMatStationaryLink(movieInfo,kalmanFilterInfoFrame1,...
%     costMatParam,nnDistFeatures,probDim,prevCost,featLifetime,...
%     trackedFeatureIndx,currentFrame)
%
%INPUT  movieInfo             : An nx1 array (n = number of frames in
%                               movie) containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features.
%                                   1st column for values and 2nd column
%                                   for standard deviations.
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%      kalmanFilterInfoFrame1 : Not relevant for this cost function.
%      costMatParam           : Structure containing variables needed for cost
%                               calculation. Contains the fields:
%             .searchRadius       : Search radius for linking features (in pixels).
%      nnDistFeatures         : Not relevant for this cost function.
%      probDim                : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%      prevCost               : Not relevant for this cost function.
%      featLifetime           : Not relevant for this cost function.
%      trackedFeatureIndx     : The matrix of feature index connectivity up
%                               to current frame.
%                               Currently not used in this cost function.
%      currentFrame           : Current frame that is being linked to the
%                               next frame.
%
%OUTPUT costMat               : Cost matrix.
%       propagationScheme     : Not relevant for this cost function.
%       kalmanFilterInfoFrame2: Not relevant for this cost function.
%       nonlinkMarker         : Value indicating that a link is not allowed.
%       errFlag               : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, August 2011

%% Output

costMat = [];
propagationScheme = [];
kalmanFilterInfoFrame2 = [];
nonlinkMarker = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatStationaryLink')
    disp('--costMatStationaryLink: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%extract the two frames of interest from movieInfo
movieInfo = movieInfo(currentFrame:currentFrame+1);

%get parameters
searchRadius = costMatParam.searchRadius;

%% Inter-particle distance calculation

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%put the coordinates of features in the 1st frame in one matrix
coord1 = movieInfo(1).allCoord(:,1:2:end);

%put the coordinates of features in the 2nd frame in one matrix
coord2 = movieInfo(2).allCoord(:,1:2:end);

%calculate the distances between features
costMat = createDistanceMatrix(coord1,coord2);

%% Search radius

%assign NaN to costs corresponding to distance > searchRadius
costMat(costMat>searchRadius) = NaN;

%square the cost matrix to make the cost = distance squared
costMat = costMat.^2;

%% Birth and death

maxCost = 2*max(max(costMat(:)),eps);

deathCost = maxCost * ones(numFeaturesFrame1,1);
birthCost = maxCost * ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
% costLR = min(min(min(costMat))-1,-1);
costLR = maxCost;
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%% nonLinkMarker

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;

%% ~~~ the end ~~~
