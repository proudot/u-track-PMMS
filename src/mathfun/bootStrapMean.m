function [conf,meanS] = bootStrapMean(variable,alpha,nBoot)
%This subFunction bootstrap the mean value of the input "variable"
% bootci calculates the confidence interval at alpha level based on the
% percentile (corrected for bias) of the distribution for the speficic
% stats
%
% Usage: [conf,meanS] = bootStrapMean(variable,alpha,nBoot)
%
% Input:
%       variable - series of points - a given stats or raw data 
%       alpha    - confidence level - see bootci
%       nBoot    - number of boostrap samples
%
% Output:
%        conf  - confidence interval for the mean value
%        meanS - mean value 
%
% Marco Vilela, 2012

opt = statset('UseParallel','never');
if matlabpool('size')
    opt = statset('UseParallel','always');
end

[conf,meanSample] = bootci(nBoot,{@nanmean,variable},'alpha',alpha,...
    'type','bca','Options',opt);

meanS = nanmean(meanSample);

end
