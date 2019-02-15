function [diffCoeffitient, pValue, mSqDisp, linSlope, yIntercept] = ...
         diffusion(positions,analyzingParameters)
     
% DIFFUSION  function computes the diffusion coeffitient, the corresponding
%            p-value, the MSD, the slope and yIntercept of the linear fit
%            for a timepoint of a given timewindow
%
% SYNOPSIS  [diffCoeffitient, pValue, mSqDisp, linSlope, yIntercept] = ...
%           diffusion(positions,analyzingParameters);
% 
% INPUT     positions =       1-by-(1-2) array of positions and covariances
%                             First column: position of tag.
%                             Second column (opt): position of reference
%               Fields:
%                   .coordinates    n-by-d array of coordinates
%                   .covariances    d-by-d-by-n array of covariances
%
%           analyzingParameters =  structure with fields:
%
%            nDiff            = number of timpoints of MSD which are used 
%                               for the linear regression fit to determine 
%                               the diffusion coeffitions
%
%            nDev             = number of timepoints of MSD compared to the
%                               linear regression fit to determine the
%                               deviation between them.
%
%            nMSD (optional)  = number of timepoints used for the 
%                               calculation
%                               of the MSD (min = 4, max = tWinSizeMin-2,
%                               default = 20) analyzingParameters   
%
% OUTPUT    diffCoeffitient  = the diffusion coeffitient of the positions
%                              calculated with the slope of the fit to the
%                              first 5 MSDs of the positions
% 
%           pValue           = the probality that the sigmaHatSquare is
%                              fischer distributed with a mean of one (a 
%                              value of 0.5 reflects 100%)
%
%           mSqDisp          = mean squared displacement of coordinates 
%                              that come with covariances
%               
%                              x-by-3 array with
%                                     - r^2=<[r(t+dt)-r(t)]^2>
%                                     - sigma(r^2) - for SEM divide by 
%                                           sqrt(nDataPoints!)
%                                     - nDataPoints
%
%           linSlope         = slope of the linear fit to the first nDev
%                              points of the MSDs
%
%           yIntercept       = the y-intercept of the linear fit to the 
%                              first nDev points of the MSDs
%
% CREATED gp 2/19/07

% calles the function meanSquaredDisplacement to compute the MSD in a given
% time window of one trajectory
mSqDisp = meanSquaredDisplacementTrackEdit(positions, [],[],analyzingParameters);

% find ~NaNs and store the location in goodData. This 
goodData = find(isfinite(mSqDisp(:,1)) & isfinite(mSqDisp(:,2)));

% creates a row vector timePoints from 1 to analyzingParameters.nDiff 
timePoints = (1:1:analyzingParameters.nDiff)';

% checks if the mean square displacement has enough values 
% (analyzingParameters.nDiff) needed for the linear fit 
if (size(goodData,1) < analyzingParameters.nDiff);
    
    % if not then the diffussion coeffitient, p-value, linSlope and 
    % yIntercept  can not be computed and their output is NaN
    diffCoeffitient = NaN;
    pValue = NaN;
    linSlope = NaN;
    yIntercept = NaN;

    % returns to the function which called this function
    return;
    
else

    %converts the output from meanSquaredDisplacementTrackEdit into the
    %variances (mSdSigma) of the MSDs
    mSdSigma = mSqDisp(1:analyzingParameters.nDiff,2)./sqrt(mSqDisp...
    (1:analyzingParameters.nDiff,3));

    % calls the function linFit which fits a linear regression to the MSD 
    % regarding its variance
    [linSlope,yIntercept,sigmaHatSquare, pValue] = linFit(timePoints,...
    mSqDisp(1:analyzingParameters.nDiff,1), mSdSigma, 2, analyzingParameters.nDiff);

    % computes the diffussion coeffitient using the slope of the linear
    % equation
    diffCoeffitient = linSlope/4;

end; % end for if else function

