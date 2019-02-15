function kalmanFilterInfo = indeKalmanResMemLMOnlineVar(numFrames,numFeatures,probDim,costMatParam)
%KALMANRESMEMLM reserves memory for Kalman filter structure for linear motion model
%
%SYNPOSIS kalmanFilterInfo = kalmanResMemLM(numFrames,numFeatures,probDim)
%
%INPUT   numFrames  : Number of frames in movie.
%        numFeatures: An array with number of feaures in each frame.
%        probDim    : Problem dimensionality.
%
%OUTPUT   kalmanFilterInfo: Structure array with number of entries equal to 
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%             .stateNoiseRSS     : RSS for online process noise estimation
%             .stateNoiseMean    : Mean for online process noise estimation
%             .stateNoiseNum     : sample number for online process noise estimation

%Philippe Roudot Dec 2012
%After Khuloud Jaqaman

%calculate vector size
vecSize = 2 * probDim;
schemeNb=mod(costMatParam.linearMotion,3)+1; 


%go over all frames and reserve memory
for iFrame = numFrames : -1 : 1

    kalmanFilterInfo(iFrame) = struct( ...
        'stateVec',zeros(numFeatures(iFrame),vecSize,schemeNb+1),...
        'stateCov',zeros(vecSize,vecSize,numFeatures(iFrame),schemeNb+1),...
        'noiseVar',zeros(vecSize,vecSize,numFeatures(iFrame),schemeNb+1),...
        'stateNoise',zeros(numFeatures(iFrame),vecSize,schemeNb+1),...
        'stateNoiseRSS',zeros(numFeatures(iFrame),2,schemeNb+1),...
        'stateNoiseMean',zeros(numFeatures(iFrame),2,schemeNb+1),...
        'stateNoiseNum',zeros(numFeatures(iFrame),2,schemeNb+1),...
        'scheme',zeros(numFeatures(iFrame),2));

end
