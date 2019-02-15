function [trackBehavior] = transBehaviorTrack...
         (trajectoryTimepointMatrix, analyzingParameters)

% TRANSBEHAVIORTRACK computes the minimum diffusion coeffitient, the 
%                    corresponding p-Value of the sigma hat square (if it
%                    is fischer distributed with a mean of 1) , the 
%                    minimum deviation between the linear fit of the MSD 
%                    and the MSD and the maximum asymmetry of all given
%                    time windows for ech time point of a given trajectory 
%
% SYNOPSIS       [trackBehavior] = transBehaviorTrack
%                (trajectoryTimepointMatrix, analyzingParameters)
% 
% INPUT   trajectoryTimepointMatrix = matrix of different timepoints (rows)
%                                     of one trajectory
%
%         analyzingParameters:        1x1 structure with fields
%           
%            tWinSizeDiffMin  = minimum size of the rolling time window
%                               used for the calculation of the diffusion
%                               coeffitient
%
%            tWinSizeDiffMax  = maximum size of the rolling time window
%                               used for the calculation of the diffusion
%                               coeffitient
%
%            tWinSizeConMin   = minimum size of the rolling time window
%                               used for the calculation of the deviation
%                               of the MSDs from their fit
%
%            tWinSizeConMax   = maximum size of the rolling time window
%                               used for the calculation of the deviation
%                               of the MSDs from their fit
%
%            tWinSizeDirMin   = minimum size of the rolling time window
%                               used for the calculation of the directed
%                               movement
%            tWinSizeDirMax   = maximum size of the rolling time window
%                               used for the calculation of the directed
%                               movement          
%
%            tWinSizeMin      = minimum window size of all min window sizes
%
%            tWinSizeMax      = maximum window size of all max window sizes
%
%            nDiff            = number of timpoints of MSD which are used 
%                               for the linear regression fit to determine 
%                               the diffusion coeffitions
%
%            nDev             = number of timepoints of MSD compared to the
%                               linear regression fit to determine the
%                               deviation between them.
%
%            nMSD (optional)  = number of timepoints used for the calculation
%                               of the MSD (min = 4, max = tWinSizeDiffMin-2,
%                               default = 20) 
%
% OUTPUT  trackBehavior   = n-by-18 matrix with
%             n rows:  timepoints of trajectories
%             columns: 1. minimum diffusion coeffitient of all different 
%                         time windows width for each time point
%                      2. p-Value, which reflects the probability that
%                         sigma Hat Square is fischer distributed
%                         with a mean of 1.
%                      3. minimum deviation of the MSD from linearity (fit 
%                         over first nDiff points of MSD) of all different 
%                         time windows width for each time point
%                      4. maximum asymmetry of all different time windows 
%                         width for each time point 
%                    5-7. empty (NaN)
%                      8. window size used for computing the min diffusion
%                         coeffitient in column 1
%                      9. window size used for computing the min deviation
%                         of the MSD from linearity in column 3
%                     10. window size used for computing the max asymmetry
%                         in column 3
%                  11-15. empty (NaN)
%                     16. actual number of data points used for the
%                         calculation of the min diffusion coeffitient
%                         (column 1)
%                     17. actual number of data points used for the
%                         calculation of the min deviation of the MSD from
%                         linearity (column 3)
%                     18. actual number of data points used for the
%                         calculation of the max asymmetry (column 4)
%
% CREATED gp 2/13/07

% determines the row and column size of the trajectory
[rowSize,columnSize] = size(trajectoryTimepointMatrix);


% deterines the number of the different time window sizes
tWinSizeNum =  ceil((analyzingParameters.tWinSizeMax-...
    analyzingParameters.tWinSizeMin+1)/2);

% creates the matrix diffCoeffitient and sets all values in the matrix to NaN
% the number of matrixes in the third dimension corresponds to the number 
% of time window sizes (tWinSizeNum)
diffCoeffitient = NaN * ones(rowSize,1,tWinSizeNum);

% creates the matrix devMsdLinReg and sets all values in the matrix to NaN
% the number of matrixes in the third dimension corresponds to the number 
% of time window sizes (tWinSizeNum)
devMsdLinReg = NaN * ones(rowSize,1,tWinSizeNum);

% creates the matrix pValue and sets all values in the matrix to NaN
% the number of matrixes in the third dimension corresponds to the number 
% of time window sizes (tWinSizeNum)
pValue = NaN * ones(rowSize,1,tWinSizeNum);

% creates the matrix asymmetry and sets all values in the matrix to NaN
% the number of matrixes in the third dimension corresponds to the number 
% of time window sizes (tWinSizeNum)
asymmetry = NaN * ones(rowSize,1,tWinSizeNum);

% creates the matrix trackBehavior and sets all values in the matrix to NaN
% the number of matrixes in the third dimension corresponds to the number 
% of time window sizes (tWinSizeNum)
trackBehavior = NaN * ones(rowSize,18);

%loops through the different time window sizes 
for i = 1:tWinSizeNum;
    
    % tWinSize increases by 2 every cycle and starts from tWinSizeMin to
    % tWinSizeMax
    tWinSize = 2*i+analyzingParameters.tWinSizeMin-2;
    
    % sets the variable winSizeMatrixCount to 1 and increases it by one
    % every cycle. It defines the third dimension of the matrix
    % trackbehavior
    winSizeMatrixCount = i;
    
    % defines the start point (row)of the rolling time window
    tWinSeriesStart = ceil(tWinSize/2);

    % defines the end point (row) of the rolling time window
    tWinSeriesEnd = rowSize-ceil(tWinSize/2)+1;

    % loops through the differemt time points (rows) in the given time 
    % window of a given trajectory
    for n = tWinSeriesStart:tWinSeriesEnd;
        
        % initialize diffYes
        diffYes = 0;

        % changes the startpoint (row)of the rolling time window during rolling
        timeWinStart = n-ceil(tWinSize/2)+1;

        % changes the endpoint (row)of the rolling time window during rolling
        timeWinEnd = timeWinStart + tWinSize - 1;

        % determines the values in the matrix trajectoryTimepointMatrix 
        % between row-possition timeWinStart and timeWinEnd for all columns
        % and stores it in the variable timeWinMatrix.
        timeWinMatrix = trajectoryTimepointMatrix(timeWinStart:timeWinEnd,:);
        
        % determines the number of data points in the timeWinMatrix
        realDataPoints(n,:,i) = (timeWinEnd-timeWinStart+1)-sum(isnan(timeWinMatrix(:,1)));
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Diffusion Coeffitient
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        % calculates the diffusion if the time window size is in between
        % tWinSizeDiffMin and twinSizeDiffMax or between tWinSizeConMin and
        % tWinSizeConMax, cause the diffCoeffitient is needed to calculate
        % the deviation which is needed to determine constrained motion
        if tWinSize >= analyzingParameters.tWinSizeDiffMin &&...
           tWinSize <= analyzingParameters.tWinSizeDiffMax;  
            %----------------------------------------------------
            % Convertion of the timeWinMatrix
            %----------------------------------------------------
            % converts the timeWinMatrix into a structure array which is the
            % input format for the meanSquaredDisplacement function.

            % coordinates include the x and y values of the timeWinMatrix 
            positions.coordinates = [timeWinMatrix(:,1),timeWinMatrix(:,2)];

            %calculates the square of sigmaX for the timeWinMatrix (column 5)
            sigmaX2 = timeWinMatrix(:,5).^2;

            %calculates the square of sigmaY for the timeWinMatrix (column 6)
            sigmaY2 = timeWinMatrix(:,6).^2;

            % finds the sizes of the matrix timeWinMatrix
            [rowSizeTimeWinMatrix, columnSizeTimeWinMatrix] = size(timeWinMatrix);  

            % creates the 2-by-2-by-rownumber stucture matrix 
            % positions.covariances and sets all values to zero. 
            positions.covariances = zeros(2,2,rowSizeTimeWinMatrix);

            % fills the structure matrix positions.covariances with the sigmaX2
            % values from the timeWinMatrix
            positions.covariances(1,1,:) = sigmaX2;

            % fills the structure matrix positions.covariances with the sigmaY2
            % values from the timeWinMatrix
            positions.covariances(2,2,:) = sigmaY2;


            % calles the function diffusion which computes the 
            % diffusion coeffitient per time point over the given window size
            [diffCoeffitient(n,:,winSizeMatrixCount),...
                pValue(n,:,winSizeMatrixCount),mSqDisp, linSlope, yIntercept]...
                = diffusion(positions,analyzingParameters);
            
            % diffusion coeffitient has been computed
            diffYes = 1;
            
        end % of diffCoeffitient calculation   
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Deviation of MSD
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                
        % calculates the deviation if the time window size is in between
        % tWinSizeConMin and twinSizeConMax 
        if tWinSize >= analyzingParameters.tWinSizeConMin &&...
                tWinSize <= analyzingParameters.tWinSizeConMax
           
            % if the diffusion coeffitient has already been computed
            if diffYes == 1;
               % calles the function standardDev which calculates the mean deviation
               % between the regression fit of the MSD (*nDiff) and the MSD of the 
               % different timepoints. 
               devMsdLinReg(n,:,winSizeMatrixCount) = linearDev...
                   (analyzingParameters,mSqDisp,linSlope, yIntercept);
               
            % if the diffusion coeffitient has not yet been computed
            else 
                
            %----------------------------------------------------
            % Convertion of the timeWinMatrix
            %----------------------------------------------------
            % converts the timeWinMatrix into a structure array which is the
            % input format for the meanSquaredDisplacement function.

            % coordinates include the x and y values of the timeWinMatrix 
            positions.coordinates = [timeWinMatrix(:,1),timeWinMatrix(:,2)];

            %calculates the square of sigmaX for the timeWinMatrix (column 5)
            sigmaX2 = timeWinMatrix(:,5).^2;

            %calculates the square of sigmaY for the timeWinMatrix (column 6)
            sigmaY2 = timeWinMatrix(:,6).^2;

            % creates the 2-by-2-by-rownumber stucture matrix 
            % positions.covariances and sets all values to zero. 
            positions.covariances = zeros(2,2,size(timeWinMatrix,1));

            % fills the structure matrix positions.covariances with the sigmaX2
            % values from the timeWinMatrix
            positions.covariances(1,1,:) = sigmaX2;

            % fills the structure matrix positions.covariances with the sigmaY2
            % values from the timeWinMatrix
            positions.covariances(2,2,:) = sigmaY2;


            % calles the function diffusion which computes the 
            % diffusion coeffitient per time point over the given window size
            [diffCoeffitientForDev, pValueForDev,mSqDisp, linSlope,...
                yIntercept] = diffusion(positions,analyzingParameters);   
           
            % calles the function standardDev which calculates the mean deviation
            % between the regression fit of the MSD (*nDiff) and the MSD of the 
            % different timepoints. 
            devMsdLinReg(n,:,winSizeMatrixCount) = linearDev...
                (analyzingParameters,mSqDisp,linSlope, yIntercept);
            
            end % ofdiffYes == 1
            
        end % of deviation calculation
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Asymmtry of track
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
        % calculates the assymetrie if the time window size is in between
        % tWinSizeDirMin and twinSizeDirMax
        if tWinSize >= analyzingParameters.tWinSizeDirMin &&...
           tWinSize <= analyzingParameters.tWinSizeDirMax
            % calles the function asymDetermination which calculates the
            % asymmetry of the gyration tensor of the x and y values of a given
            % time window
            asymmetry(n,:,winSizeMatrixCount) = asymDetermination(timeWinMatrix);
        end % of asymmetry calculation
        
    end % of time point loop
    
    % to check the progress of the computation during running
    %disp(winSizeMatrixCount);
    
end % of time window size loop

% finds the minimum diffusion coeffitients of all time window sizes for
% each timepoint which is stored in column 1 of the matrix trackBehavior.
% Also creates the Vector minDiffCoeffitientIndex which has the Index
% information in which third dimension of the diffCoeffitient
% the minimum is found. This information is used to find the corresponding
% p-Values.
[trackBehavior(:,1), trackBehavior(:,8)] = min(diffCoeffitient,[],3);

% loops through all minDiffCoeffitiens and adds the corresponding p-value
% to each minDiffCoeffitient-value in the second column of the matrix 
%trackBehavior corresponding to the minDiffCoeffitientIndex
for w = tWinSeriesStart-i+1:tWinSeriesEnd+i-1
    trackBehavior(w,2) = pValue(w,1,trackBehavior(w,8));
    trackBehavior(w,16) = realDataPoints(w,1,trackBehavior(w,8));
end

% finds the minimum deviations between the regression fit of the MSD 
% (*nDiff) and the MSD of the different timepoints of all time window sizes
% which is stored in column 3 of the matrix trackBehavior and saves the
% corresponding window size in in column 7 of the matrix trackBehavior for
% the later calculation of the duration of the transient behavior. 
[trackBehavior(:,3),trackBehavior(:,9)] = min(devMsdLinReg,[],3);

% loops through all minimum deviations and adds the corresponding real data
% points to each minDevMsdLinReg-value in the 14 column of the matrix 
% trackBehavior.
for w = tWinSeriesStart-i+1:tWinSeriesEnd+i-1
    trackBehavior(w,17) = realDataPoints(w,1,trackBehavior(w,9));
end

% finds the maximum asymmetry of the x and y values (possitions) of the 
% different timepoints of all time window sizes which is stored in the 
% column 4 of the matrix trackBehavior and saves the
% corresponding window size in column 8 of the matrix trackBehavior for the
% later calculation of the duration of the transient behavior.  
[trackBehavior(:,4), trackBehavior(:,10)] = max(asymmetry,[],3);

% loops through all maximum asymmetries and adds the corresponding real data
% points to each maximum asymmetry-values in the 15 column of the matrix 
% trackBehavior.
for w = tWinSeriesStart-i+1:tWinSeriesEnd+i-1
    trackBehavior(w,18) = realDataPoints(w,1,trackBehavior(w,10));
end

% correction of the index of the parameters minDiff, minDev and maxAsym  
trackBehavior(:,8:10) = (trackBehavior(:,8:10)).*2-2+ analyzingParameters.tWinSizeMin;