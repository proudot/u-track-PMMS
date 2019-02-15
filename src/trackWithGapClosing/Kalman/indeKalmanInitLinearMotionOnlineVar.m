function [indeKalmanFilterInfo,errFlag] = indeKalmanInitLinearMotionOnlineVar(frameInfo,...
    probDim,costMatParam)
%KALMANINITLINEARMOTION initializes Kalman filter state vector and covariance matrix for features in a frame
%
%SYNOPSIS [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo,...
%    probDim,costMatParam)
%
%INPUT  frameInfo       : Structure with fields (for 1 frame):
%             .allCoord     : Image coordinate system x,dx,y,dy,[z,dz] of 
%                             detected features and their uncertainties (in
%                             pixels).
%             .num          : Number of features in frame.
%       probDim         : Problem dimension. 2 (for 2D) or 3 (for 3D).
%       costMatParam    : Linking cost matrix parameters. In addition
%                         to the fields needed for costMatLinearMotionLink,
%                         if it has the field 
%             .kalmanInitParam, then sub-fields within will be used for
%                               initialization. 
%                   .convergePoint: Convergence point (x, y, [z]
%                                   coordinates) of tracks if motion is
%                                   radial, in image coordinate system.
%                                   Should be a row vector. If supplied,
%                                   radial form is assumed and value is
%                                   used to estimate initial velocities. If
%                                   not supplied, then initial velocity is
%                                   taken as zero.
%                   .searchRadiusFirstIteration: max number of pixels to
%                                   use in the first frame pair where the
%                                   search is essentially a nearest
%                                   neighbor. here we use this parameter to
%                                   get the initial noise variance. in
%                                   kathryn's alternate function
%                                   costMatLinearMotionLink_EB3.m, the
%                                   min/max search radii are not applied in
%                                   the first iteration so that this
%                                   initial radius can be sufficiently big
%                                   if features are moving quickly (ie
%                                   faster than the max search radius)
%                   .initVelocity : Initial velocity guess (vx, vy, [vz]),
%                                   in whatever units coordinates are in
%                                   per frame.
%
%OUTPUT kalmanFilterInfo: Structure with fields:
%             .stateVec     : State vector for each feature.
%             .stateCov     : State covariance matrix for each feature.
%             .noiseVar     : Variance of state noise for each feature.
%             .stateNoiseRSS     : RSS for online process noise estimation
%             .stateNoiseMean    : Mean for online process noise estimation
%             .stateNoiseNum     : sample number for online process noise estimation
%             .scheme:
%       errFlag         : 0 if function executes normally, 1 otherwise.
%
%
%Philippe Roudot Dec 2012 
%after Khuloud Jaqaman, March 2007

[kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo, probDim,costMatParam);

%find number of features in frame
numFeatures = frameInfo.num;


indeKalmanFilterInfo = kalmanFilterInfo;

numScheme= mod(costMatParam.linearMotion,3)+1;

% Initialize state vector for each indepenant KF (+1 for the general
% estimate)
for i = 1:numScheme
indeKalmanFilterInfo.stateVec = cat(3,indeKalmanFilterInfo.stateVec,kalmanFilterInfo.stateVec);
indeKalmanFilterInfo.stateCov = cat(4,indeKalmanFilterInfo.stateCov,kalmanFilterInfo.stateCov);
indeKalmanFilterInfo.noiseVar = cat(4,indeKalmanFilterInfo.noiseVar,kalmanFilterInfo.noiseVar);
end

%Initialize online process noise estimation variable
indeKalmanFilterInfo.stateNoiseRSS=zeros(size(kalmanFilterInfo.stateVec,1),2,numScheme+1);
indeKalmanFilterInfo.stateNoiseMean=zeros(size(kalmanFilterInfo.stateVec,1),2,numScheme+1);
indeKalmanFilterInfo.stateNoiseNum=zeros(size(kalmanFilterInfo.stateVec,1),2,numScheme+1);



%% %%%% ~~ the end ~~ %%%%
