function cumDistrExpo = calcCumDistrExpo(param,abscissa)
%CALCCUMDISTREXPO calculates the cumulative distribution of an exponential
%
%SYNOPSIS cumDistrExpo = calcCumDistrExpo(param,abscissa)
%
%INPUT  param         : Vector of parameters indicating the exponent and
%                       amplitude of the exponential.
%       abscissa      : Abscissa values at which the cumulative
%                       distribution is calculated.
%
%OUTPUT cumDistrExpo  : Values of the resulting cumulative distribution
%                       given the input abscissa values.

%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cumDistrExpo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--calcCumDistrExpo: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the cumulative distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the cumulative distribution
cumDistrExpo = param(2)*expcdf(abscissa,param(1));

%%%%% ~~ the end ~~ %%%%%

