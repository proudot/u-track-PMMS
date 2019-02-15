function [numObsPerBin,binCenter,expoParam,errFlag] = fitHistWithExponential(...
    observations,showPlot)
%FITHISTWITHEXPONENTIAL fits an exponential to a histogram
%
%SYNOPSIS [numObsPerBin,binCenter,expoParam,errFlag] = fitHistWithExponential(...
%    observations,showPlot)
%
%INPUT  observations: Vector of observations whose histogram is to be fitted.
%       showPlot    : 1 to plot the histogram and fitted exponential, 0 otherwise. 
%                     Optional. Default: 1.
%
%OUTPUT numObsPerBin: Number of observations that fall in each bin.
%       binCenter   : Center of each bin.
%       expoParam   : Row vector where 1st entry is the exponent and 2nd 
%                     entry is the amplitude.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS The fitted exponential is normalized. Thus, the fitted equation is
%given by (expoParam(2)/expoParam(1))*exp(-x/expoParam(1))
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numObsPerBin = [];
binCenter = [];
expoParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--fitHistWithExponential: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(showPlot)
    showPlot = 1;
else
    if showPlot ~= 0 && showPlot ~= 1
        disp('fitHistWithExponential: Variable "showPlot" should be either 0 and 1!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Histogram calculation and fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the number of observations
numObservations = length(find(~isnan(observations)));

%calculate the histogram
[numObsPerBin,binCenter] = optimalHistogram(observations);
numObsPerBin = numObsPerBin';
binCenter = binCenter';

%determine the number of bins used
numBins = length(binCenter);

%calculate the cumulative histogram
cumHist = zeros(numBins,1);
for iBin = 1 : numBins
    cumHist(iBin) = sum(numObsPerBin(1:iBin));
end

%set some optimization options
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-6);

%assign initial values to unknown parameters
x0 = [1 numObservations/2]';

%assign lower bounds
lb = zeros(2,1);

%assign upper bounds
ub = [];

%estimate unknown parameters
expoParam = lsqcurvefit(@calcCumDistrExpo,x0,binCenter,cumHist,lb,ub,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the user wants to plot
if showPlot

    %get the distribution from the optimized parameters
    distrExpo = expoParam(2)*exppdf(binCenter,expoParam(1))*(binCenter(end)-binCenter(end-1));

    %get the cumulative distribution from the optimized parameters
    cumDistrExpo = expoParam(2)*expcdf(binCenter,expoParam(1));

    %make new figure
    figure

    %plot the histogram and the fitted exponential in the left half of the
    %figure
    subplot(1,2,1);
    plot(binCenter,numObsPerBin,'k.')
    hold on
    plot(binCenter,distrExpo,'r')

    %plot the histogram and the fitted exponential in the right half of the
    %figure
    subplot(1,2,2);
    plot(binCenter,cumHist,'k.')
    hold on
    plot(binCenter,cumDistrExpo,'r')

end %(if showPlot)

%%%%% ~~ the end ~~ %%%%%
