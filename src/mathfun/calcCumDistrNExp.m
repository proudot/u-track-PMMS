function cumDistrNExp = calcCumDistrNExp(param,abscissa)
%CALCCUMDISTRNEXP calculates the cumulative distribution of N exponentials
%
%SYNOPSIS cumDistrNExp = calcCumDistrNExp(param,abscissa)
%
%INPUT  param         : Vector of parameters indicating the means and
%                       amplitude of the exponentials.
%       abscissa      : Abscissa values at which the cumulative
%                       distribution is calculated.
%
%OUTPUT cumDistrNExp  : Values of the resulting cumulative distribution
%                       given the input abscissa values.
%       jacobianMat   : Jacobian matrix.

%Khuloud Jaqaman, June 2012

%% Output

cumDistrNExp = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--calcCumDistrNExp: Incorrect number of input arguments!');
    return
end

%% Calculating the cumulative distribution

%get number of exponentials
numExp = length(param)/2;

%get their means and amplitudes
expMean = param(1:numExp);
expAmp  = param(numExp+1:end);

%get number of data points
numData = length(abscissa);

%calculate the cumulative distribution and the jacobian matrix
cumDistrNExp = zeros(numData,1);
% jacobianMat = zeros(numData,2*numExp);
for iExp = 1 : numExp
    cumDistrNExp = cumDistrNExp + expAmp(iExp)*(1 - exp(-abscissa/expMean(iExp)));
    %     jacobianMat(:,iExp) = -(expAmp(iExp)/expMean(iExp)^2)*exp(-abscissa/expMean(iExp));
    %     jacobianMat(:,numExp+iExp) = 1 - exp(-abscissa/expMean(iExp));
end
    
%% ~~ the end ~~
    
