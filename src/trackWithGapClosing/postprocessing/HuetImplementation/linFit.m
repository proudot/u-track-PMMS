function [a,b,sigmaHatSquare, pValue] = linFit(x, y, mSdSigma, u, nData)

% LINFIT fits a linear regression to datapoints(x,y) regarding their
%        variance (mSdSigma) and computes the probability (pValue) that the
%        sigma hat squares (sigmaHatSquare) are fischer distributed with a 
%        mean of one
%
% SYNOPSIS [a,b,sigmaHatSquare, pValue] = linFit(x, y, sigma, u, nData)
%
% INPUT    x        = x-values (data) for the linear fit
%                  
%          y        = y-values (data) for the linear fit  
%
%          mSdSigma = variance of the data (MSD)
% 
%          u        = number of parameters
%
%          nData    = number of datapoints for fischer distribution
%
% OUTPUT   a             = slope of the linear fit (linSlope)
%
%          b             = the y-intercept of the linear fit (yIntercept)
%
%          sigmaHatSuare = sigma hat square, which reflects the minimum
%                          average deviation between the fitted line to the
%                          data points and each data point taking the  
%                          variance of each data point in account
%         
%          pValue        = the probality that the sigmaHatSquare is
%                          fischer distributed with a mean of one (a value
%                          of 0.5 reflects 100%)
%
% CREATED  gp 03/02/07

% determines the row and column size of the vector with the x values
[rowSizeX, columnSizeX] = size(x);

% defining the Sums for the calculation of the slope (linSlope) and the
% y-intercept (yIntercept)
S = sum(1./mSdSigma.^2);
Sx = sum(x./mSdSigma.^2);
Sy = sum(y./mSdSigma.^2);
Sxx = sum((x.^2)./(mSdSigma.^2));
Sxy = sum((x.*y)./mSdSigma.^2);

% defines the variable (delta) which is used for the simplification of the
% calculation of linSlope and yIntercept
delta = S*Sxx - (Sx)^2;

% computes the slope (linSlope) of the linear regression of the MSDs 
a = (S*Sxy - Sx*Sy)/delta;

% computes the Y-intercept (yIntercept)of the linear regression of the MSDs
b = (Sxx*Sy - Sx*Sxy)/delta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additional computations which were in the original paper but are not
% necessary for the usage of this function
%
% % computes the varience (sigmaLinSlope) of the slope (linSlope) of the
% % linear regression of the MSD
% sigmaLinSlope = sqrt(S/delta);
% 
% % computes the varience (sigmaYIntercept) of the Y-intercept (yIntercept)
% % of the linear regression of the MSD
% sigmaYIntercept = sqrt(Sxx/delta);
% 
% % computes the covariance (covLinSlopeInter) of (linSlope) and (yIntercept)
% covLinSlopeInter = (-Sx)/delta;
% 
% % computes the correlation coeffitient (rLinSLopeYIntercept) between the
% % uncertainty in (yIntercept) and the uncertainty in (LinSlope) which is a 
% % number between -1 and 1. A positive value of (rLinSlopeYIntercept) 
% % indicates that the errors in both variables have the same sign. A
% % negative value indicates that the errors are anticorrelated and are
% % likely to have opposite signs.
% rLinSlopeYIntercept = (-Sx)/(sqrt(S*Sxx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creates the probability matrix P, which contains the 1/variance 
% (diagonal)of the data points (MSD). The 1 is used cause we assume that 
% the sigmaHatSquares are fischer distributed with a mean of 1.
P = diag(1./mSdSigma.^2);

% calculates the y-values (yFunc) of a linear distribution with the
% computed a (linSlope) and b (yIntercept) values and the x-values from the
% data (MSD). This yFunc represents the model of the fitted data.
yFunc = a * x + b;

% computes the difference between the original data points (y-values) and
% the computed y-values (yFunc). Difference between the real data and the
% model.
e = y-yFunc;

% chiSquare = e'*P*e (sum of differences between the model (calculated with the
% computed a and b value) and the meassurements weighted by their variance.
% sigmaHatSquare is chiSquare divided by the degrees of freedom. (min mean
% of the difference of each meassurement to the model in y)
sigmaHatSquare = (e'*P*e)/(rowSizeX-u);

% computes the goodness-of-fit of the MSD-data to the model. To verify the
% goodness of the fit, it is assumed that sigmaHatSquare is distributed by
% a fischer distribution with a mean of 1. The pValue reflects the
% probability that the caculated sigmaHatSquare belongs to this
% distribution. A pValue of 0.5 reflects a 100% probability that the
% sigmaHatsquare is distributed by fischer with mean 1. 
if sigmaHatSquare <= 1
    pValue = fcdf(sigmaHatSquare,rowSizeX-u,Inf);
else
    pValue = 1-(fcdf(sigmaHatSquare,rowSizeX-u,Inf));
end




  
