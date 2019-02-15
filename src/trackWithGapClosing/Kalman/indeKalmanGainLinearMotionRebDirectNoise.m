function [kalmanFilterInfoOut,errFlag] = indeKalmanGainLinearMotionRebDirectNoise(trackedFeatureIndx,...
    frameInfo,kalmanFilterInfoTmp,propagationScheme,kalmanFilterInfoIn,probDim,...
    filterInfoPrev,costMatParam,initFunctionName)
%KALMANGAINLINEARMOTION uses the Kalman gain from linking to get better estimates of the state vector, its covariance matrix, state noise and its variance
%
%SYNOPSIS [kalmanFilterInfoOut,errFlag] = kalmanGainLinearMotion(trackedFeatureIndx,...
%    frameInfo,kalmanFilterInfoTmp,propagationScheme,kalmanFilterInfoIn,probDim,...
%    filterInfoPrev,costMatParam,initFunctionName)
%
%INPUT  trackedFeatureIndx : Matrix showing connectivity between features
%                            in current frame (listed in last column of matrix) 
%                            and features in previous frames. A zero in 
%                            columns before last indicates that feature is
%                            not connected to any previous features.
%       frameInfo          : Structure with fields (for current frame):
%             .allCoord        : x,dx,y,dy,[z,dz] of features collected in one
%                                matrix.
%             .amp             : Amplitudes of PSFs fitting detected features. 
%                                1st column for values and 2nd column 
%                                for standard deviations.
%       kalmanFilterInfoTmp: Structure with fields (for current frame):
%             .stateVec        : Kalman filter prediction of state
%                                vector of all features in current frame 
%                                based on all 3 motion models.
%             .stateCov        : Kalman filter prediction of state
%                                covariance matrix of all features in
%                                current frame based on all 3 motion models.
%             .obsVec          : Kalman filter prediction of the
%                                observed variables for all features in 
%                                current frame based on all 3 motion models.
%       propagationScheme  : Matrix indicating the propagation scheme that
%                            yielded the lowest cost for a link between two
%                            features.
%       kalmanFilterInfoIn : Structure with fields (for all previous frames):
%             .stateVec        : Kalman filter state vector for all features.
%             .stateCov        : Kalman filter state covariance matrix
%                                for all features.
%             .noiseVar        : Variance of state noise for all feature.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       probDim            : Problem dimension. 2 (for 2D) or 3 (for 3D).
%       filterInfoPrev     : Structure with fields (for current frame):
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%                            Enter [] if there is no previous information.
%       costMatParam       : Linking cost matrix parameters.
%       initFunctionName   : Name of function for Kalman filter
%                            initialization.
%
%OUTPUT kalmanFilterInfoOut: Structure with fields (for all frames upto current):
%             .stateVec        : Kalman filter state vector for all features.
%             .stateCov        : Kalman filter state covariance matrix
%                                for all features.
%             .noiseVar        : Variance of state noise for all features.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       errFlag            : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2007

%% Output

kalmanFilterInfoOut = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 9
    disp('--kalmanGainLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

if isempty(filterInfoPrev)
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%% Gain calculation and update

%take absolute value of all noise variances - this takes care of the
%negative variances used to indicate first appearances
for iFrame = 1 : length(kalmanFilterInfoIn)
    kalmanFilterInfoIn(iFrame).noiseVar = abs(kalmanFilterInfoIn(iFrame).noiseVar);
end

%copy kalmanFilterInfoIn into kalmanFilterInfoOut
kalmanFilterInfoOut = kalmanFilterInfoIn;

%get number of features in current frame and frame number
[numFeatures,iFrame] = size(trackedFeatureIndx);

%construct Kalman filter observation matrix
observationMat = [diag(ones(probDim,1)) zeros(probDim)];

%calculate vector sizes
vecSize = 2 * probDim;


switch costMatParam.linearMotion
    
    case 0 %only random motion
        
        transMat(:,:,1) = eye(vecSize); %zero drift transition matrix
        numSchemes = 1;
        
    case 1 %random motion + directed motion
        
        transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
        transMat(:,:,2) = eye(vecSize); %zero drift transition matrix
        numSchemes = 2;
        
    case 2 %random motion + directed motion that can switch to opposite direction at any moment
        
        transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
        transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
        transMat(:,:,3) = eye(vecSize); %zero drift transition matrix
        numSchemes = 3;
    case 3 %random motion + directed motion that can switch to opposite direction at any moment
        
        transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
        numSchemes = 1;
        
end

%construct observation matrix
observationMat = [eye(probDim) zeros(probDim)]; %observation matrix


%go over all features in current frame
for iFeature = 1 : numFeatures

    %find index of feature in previous frame that this feature is connected to
    iFeaturePrev = trackedFeatureIndx(iFeature,end-1);

    %if this feature is connected to a feature in previous frame
    if iFeaturePrev ~= 0
        % A "main" KF is updated to estimate process noise variance involved in cutoff estimation.
        
        % find propagation scheme leading to this link and save in kalmanFilterInfo
        bestScheme = propagationScheme(iFeaturePrev,iFeature);
        kalmanFilterInfoOut(iFrame).scheme(iFeature,1) = bestScheme; %to current feature
        kalmanFilterInfoOut(iFrame-1).scheme(iFeaturePrev,2) = bestScheme; %from previous feature
                
        %get the best prediction and observation vector
        stateVecOld = kalmanFilterInfoIn(iFrame-1).stateVec(iFeaturePrev,:,numSchemes+1)' ;
        obsVecOld = observationMat*stateVecOld;
        
        %Use the main KF covariance
        stateCovOld = kalmanFilterInfoOut(iFrame-1).stateCov(:,:,iFeaturePrev,numSchemes+1);
        noiseVar = kalmanFilterInfoOut(iFrame-1).noiseVar(:,:,iFeaturePrev,numSchemes+1);

        %predict state covariance matrix of feature in 2nd frame
        stateCovPred = transMat(:,:,numSchemes- mod(bestScheme,  numSchemes)) ... 
            *stateCovOld*transMat(:,:,numSchemes- mod(bestScheme,  numSchemes))' ...
            + noiseVar;
        
        %calculate Kalman gain
        kalmanGain = stateCovPred*observationMat' / ...
            (observationMat*stateCovPred*observationMat'+...
            diag(eps+frameInfo.allCoord(iFeature,2:2:end).^2)); %add epsilon to avoid division by zero in the extreme case of all errors = zero
        
        %estimate state noise in previous frame and save in kalmanFilterInfo
        stateNoise = kalmanGain * (frameInfo.allCoord(iFeature,1:2:end)' - obsVecOld);
        kalmanFilterInfoOut(iFrame-1).stateNoise(iFeaturePrev,:,numSchemes+1) = stateNoise';
        
        %update estimate of state vector in current frame
        stateVec = stateVecOld + stateNoise;
        
        %update estimate of state covariance matrix in current frame
        stateCov = stateCovPred - kalmanGain*observationMat*stateCovPred;
        
        % online variance estimation initialization
        if(bestScheme<=numSchemes)
            if((bestScheme==numSchemes)&& ... % On transition to brownian motion
               (kalmanFilterInfoOut(iFrame-1).scheme(iFeaturePrev,1)~=numSchemes)) 
                % Reboot the brownian scheme process noise. 
                stateNoiseRSS=kalmanFilterInfoIn(iFrame-1).stateNoiseRSS(iFeaturePrev,:,numSchemes);
                stateNoiseMean=kalmanFilterInfoIn(iFrame-1).stateNoiseMean(iFeaturePrev,:,numSchemes);
                stateNoiseNum=kalmanFilterInfoIn(iFrame-1).stateNoiseNum(iFeaturePrev,:,numSchemes);            
            else 
                % otherwise 
                stateNoiseRSS=kalmanFilterInfoIn(iFrame-1).stateNoiseRSS(iFeaturePrev,:,numSchemes+1);
                stateNoiseMean=kalmanFilterInfoIn(iFrame-1).stateNoiseMean(iFeaturePrev,:,numSchemes+1);
                stateNoiseNum=kalmanFilterInfoIn(iFrame-1).stateNoiseNum(iFeaturePrev,:,numSchemes+1);
            end 
        else
            % if best prediction stems from the last tracking round use its corresponding 
            % noise process online estimation variable.
            stateNoiseRSS = filterInfoPrev.stateNoiseRSS(iFeature,:,numSchemes+1);
            stateNoiseMean = filterInfoPrev.stateNoiseMean(iFeature,:,numSchemes+1);
            stateNoiseNum = filterInfoPrev.stateNoiseNum(iFeature,:,numSchemes+1);
        end 

        % Online variance estimation.
        stateNoiseMeanNp1=stateNoiseMean+(stateNoise(1:probDim:end)'- stateNoiseMean)./ (stateNoiseNum+1); 
        stateNoiseMeanNp2=stateNoiseMeanNp1+(stateNoise(2:probDim:end)'- stateNoiseMeanNp1)./(stateNoiseNum+2); 
        stateNoiseRSS=stateNoiseRSS+(stateNoise(1:probDim:end)'-stateNoiseMeanNp1).*(stateNoise(1:probDim:end)'-stateNoiseMean) + ...
            (stateNoise(2:probDim:end)'-stateNoiseMeanNp2).*(stateNoise(2:probDim:end)'-stateNoiseMeanNp1);
        stateNoiseNum=stateNoiseNum+2;
        stateNoiseMean=stateNoiseMeanNp2;
        noiseVar = zeros(1,2*probDim);
        noiseVar(1:probDim) = stateNoiseRSS(1)/(stateNoiseNum(1)-1); 
        noiseVar(probDim+1:2*probDim) = stateNoiseRSS(2)/(stateNoiseNum(2)-1);
        noiseVar=diag(noiseVar);
        
        %save this information in kalmanFilterInfo
        kalmanFilterInfoOut(iFrame).stateVec(iFeature,:,numSchemes+1) = stateVec';
        kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature,numSchemes+1) = stateCov;
        kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature,numSchemes+1) = noiseVar;
        kalmanFilterInfoOut(iFrame).stateNoiseRSS(iFeature,:,numSchemes+1) = stateNoiseRSS;
        kalmanFilterInfoOut(iFrame).stateNoiseMean(iFeature,:,numSchemes+1) = stateNoiseMean;
        kalmanFilterInfoOut(iFrame).stateNoiseNum(iFeature,:,numSchemes+1) = stateNoiseNum;
        
        % For each scheme an independant KF is updated
        for iScheme = 1 : numSchemes
            
            
            if (iScheme~=numSchemes)&&(bestScheme>numSchemes)&&(iScheme==(numSchemes- mod(bestScheme,  numSchemes))) 
                % If we used a past prediction then the past gain and other KF variable should
                % replace the current variable associated to scheme <iScheme>
                stateVecOld = kalmanFilterInfoTmp.stateVec(iFeaturePrev,:,bestScheme)';
                stateCovOld = kalmanFilterInfoTmp.stateCov(:,:,iFeaturePrev,bestScheme);
                obsVecOld = kalmanFilterInfoTmp.obsVec(iFeaturePrev,:,bestScheme)';
                
                stateNoiseRSS = filterInfoPrev.stateNoiseRSS(iFeature,:,iScheme);
                stateNoiseMean = filterInfoPrev.stateNoiseMean(iFeature,:,iScheme);
                stateNoiseNum = filterInfoPrev.stateNoiseNum(iFeature,:,iScheme);
                noiseVar=filterInfoPrev.noiseVar(:,:,iFeature,iScheme);
            else 
                stateVecOld = kalmanFilterInfoTmp.stateVec(iFeaturePrev,:,iScheme)';
                stateCovOld = kalmanFilterInfoTmp.stateCov(:,:,iFeaturePrev,iScheme);
                obsVecOld = kalmanFilterInfoTmp.obsVec(iFeaturePrev,:,iScheme)';
                
                stateNoiseRSS=kalmanFilterInfoIn(iFrame-1).stateNoiseRSS(iFeaturePrev,:,iScheme);
                stateNoiseMean=kalmanFilterInfoIn(iFrame-1).stateNoiseMean(iFeaturePrev,:,iScheme);
                stateNoiseNum=kalmanFilterInfoIn(iFrame-1).stateNoiseNum(iFeaturePrev,:,iScheme);
                noiseVar=kalmanFilterInfoIn(iFrame-1).noiseVar(:,:,iFeaturePrev,iScheme);
            end 

            %calculate Kalman gain
            kalmanGain = stateCovOld*observationMat' / ...
                (observationMat*stateCovOld*observationMat'+...
                 diag(eps+frameInfo.allCoord(iFeature,2:2:end).^2)); ...
            %add epsilon to avoid division by zero in the extreme case of all errors = zero
        
            %estimate state noise in previous frame and save in kalmanFilterInfo
            stateNoise = kalmanGain * (frameInfo.allCoord(iFeature,1:2:end)' - obsVecOld);
            kalmanFilterInfoOut(iFrame-1).stateNoise(iFeaturePrev,:,iScheme) = stateNoise';
        
            %update estimate of state vector in current frame
            stateVec = stateVecOld + stateNoise;
        
            %update estimate of state covariance matrix in current frame
            stateCov = stateCovOld - kalmanGain*observationMat*stateCovOld;

            % When this scheme is selected the process noise variance is updated.                 
            if((iScheme==bestScheme))
                % Online mean estimation of \bar x_{n+2} and \bar x_{n+1} for posistion and velocity.
                stateNoiseMeanNp1=stateNoiseMean+(stateNoise(1:probDim:end)'- stateNoiseMean)./ (stateNoiseNum+1); 
                stateNoiseMeanNp2=stateNoiseMeanNp1+(stateNoise(2:probDim:end)'- stateNoiseMeanNp1)./(stateNoiseNum+2); 
                stateNoiseRSS=stateNoiseRSS+(stateNoise(1:probDim:end)'-stateNoiseMeanNp1).*(stateNoise(1:probDim:end)'-stateNoiseMean) + ...
                    (stateNoise(2:probDim:end)'-stateNoiseMeanNp2).*(stateNoise(2:probDim:end)'-stateNoiseMeanNp1);
                stateNoiseNum=stateNoiseNum+2;                
                
                stateNoiseMean=stateNoiseMeanNp2;

                noiseVar = zeros(1,2*probDim);
                noiseVar(1:probDim) = stateNoiseRSS(1)/(stateNoiseNum(1)-1);
                noiseVar(probDim+1:2*probDim) = stateNoiseRSS(2)/(stateNoiseNum(2)-1);
                noiseVar=diag(noiseVar);
            end

            %save this information in kalmanFilterInfo
            kalmanFilterInfoOut(iFrame).stateVec(iFeature,:,iScheme) = stateVec';
            kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature,iScheme) = stateCov;
            kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature,iScheme) = noiseVar;
            kalmanFilterInfoOut(iFrame).stateNoiseRSS(iFeature,:,iScheme) = stateNoiseRSS;
            kalmanFilterInfoOut(iFrame).stateNoiseMean(iFeature,:,iScheme) = stateNoiseMean;
            kalmanFilterInfoOut(iFrame).stateNoiseNum(iFeature,:,iScheme) = stateNoiseNum;   
        end
        
       
    else %if this feature is not connected to anything in previous frame

        %initialize Kalman filter for this feature
        if usePriorInfo %use prior information if supplied
            for iScheme=1:(numSchemes+1)
                kalmanFilterInfoOut(iFrame).stateVec(iFeature,:,iScheme) = filterInfoPrev.stateVec(iFeature,:,iScheme);
                kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature,iScheme) = filterInfoPrev.stateCov(:,:,iFeature,iScheme);
                kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature,iScheme) = filterInfoPrev.noiseVar(:,:,iFeature,iScheme);
                kalmanFilterInfoOut(iFrame).stateNoiseRSS(iFeature,:,iScheme) = filterInfoPrev.stateNoiseRSS(iFeature,:,iScheme);
                kalmanFilterInfoOut(iFrame).stateNoiseMean(iFeature,:,iScheme) = filterInfoPrev.stateNoiseMean(iFeature,:,iScheme);
                kalmanFilterInfoOut(iFrame).stateNoiseNum(iFeature,:,iScheme) = filterInfoPrev.stateNoiseNum(iFeature,:,iScheme);
            end
        else
            featureInfo.allCoord = frameInfo.allCoord(iFeature,:);
            featureInfo.num = 1;
            eval(['[filterTmp,errFlag] = ' initFunctionName '(featureInfo,probDim,'...
                'costMatParam);']);
            kalmanFilterInfoOut(iFrame).stateVec(iFeature,:,:) = filterTmp.stateVec;
            kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature,:) = filterTmp.stateCov;
            kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature,:) = filterTmp.noiseVar;
            kalmanFilterInfoOut(iFrame).stateNoiseRSS(iFeature,:,:) = filterTmp.stateNoiseRSS;
            kalmanFilterInfoOut(iFrame).stateNoiseMean(iFeature,:,:) = filterTmp.stateNoiseMean;
            kalmanFilterInfoOut(iFrame).stateNoiseNum(iFeature,:,:) = filterTmp.stateNoiseNum;
            
        end

    end %(if iFeaturePrev ~= 0 ... else ...)

    %under all circumstances, fill in the real position of each particle
    for iScheme = 1 : (numSchemes+1)
        tmpState = kalmanFilterInfoOut(iFrame).stateVec(iFeature,:,iScheme);
        tmpState(1:probDim) = frameInfo.allCoord(iFeature,1:2:end);
        kalmanFilterInfoOut(iFrame).stateVec(iFeature,:,iScheme) = tmpState;
    end
    
end %(for iFeature = 1 : numFeatures)


%% %%%% ~~ the end ~~ %%%%
