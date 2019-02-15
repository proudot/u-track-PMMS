function [devMsdLinReg] = linearDev(analyzingParameters,mSqDisp,linSlope,...
                                      yIntercept)

% STANDARDDEV function that computes the standard deviation between MSD and
%             the linear regression
%
% SYNOPSIS [devMsdLinReg] = linearDev(analyzingParameters,mSqDisp,
%                                       linSlope, yIntercept)  
% 
% INPUT          analyzingParameters = structure including the fields
%                                      tWinSizeMin, tWinSizeMax, nDiff and
%                                      nDev
%
%                mSqDisp             = mean squared displacement of coordinates 
%                                      that come with covariances
%               
%                                      x-by-3 array with
%                                         - r^2=<[r(t+dt)-r(t)]^2>
%                                         - sigma(r^2) - for SEM divide by 
%                                           sqrt(nDataPoints!)
%                                         - nDataPoints
%
%                linSlope            = slope of the linear fit to the MSD
%
%                yIntercept          = Y-intercept of the linear fit to the
%                                      MSD
%
% OUTPUT         devMsdLinReg        = deviation between MSD and the linear
%                                      regression
%
% CREATED gp 2/20/07


% creates a row vector timepoints from 1 to analyzingParameters.nDiff 
timepoints = (1:1:analyzingParameters.nDev)';

% computes the y-values to the timepoints using the linear equation from
% the polynomal fit

msdDiffusion = linSlope * timepoints + yIntercept;

% defines the row and column sizes of the mSqDisp matrix
[rowSizeMSqDisp, columnSizeMSqDisp] = size(mSqDisp);

% checks if the mean square displacement has enough values 
% (analyzingParameters.nDev) needed for the computation of the deviation 
if (rowSizeMSqDisp < analyzingParameters.nDev);
    
    % if not the deviation is set to NaN
    devMsdLinReg = NaN;
    
    % returns to the function which called this function
    return;
    
else

% calculates the average deviation between the computet MSDs and the linear
% fit of the first points (nDiff) of the MSD
devMsdLinReg = mean(((mSqDisp(1:analyzingParameters.nDev,1))-msdDiffusion)./msdDiffusion);

% end for if else function
end