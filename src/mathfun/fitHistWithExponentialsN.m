function [expParam,errFlag] = fitHistWithExponentialsN(observations,alpha,...
    showPlot,maxNumExp,binStrategy,numBinIn,plotName,meanLB)
%FITHISTWITHEXPONENTIALSN determines the number of exponentials + their characteristics to fit a histogram
%
%SYNOPSIS [numObsPerBinP,binCenterP,expParam,errFlag] = fitHistWithExponentialsN(...
%    observations,alpha,showPlot,maxNumExp,binStrategy,plotName,meanLB)
%
%INPUT  observations: Vector of observations whose histogram is to be fitted.
%       alpha       : Alpha-value for the statistical test that compares the
%                     fit of n+1 exponentials to the fit of n exponentials.
%                     Optional. Default: 0.01.
%       showPlot    : 0 to not plot anything.
%                     1 to plot the histogram and fitted exponentials.
%                     2 as 1, but with smooth histogram.
%                     Optional. Default: 1.
%       maxNumExp   : Upper limit on the number of exponentials.
%                     Optional. Default: 100.            
%       binStrategy : Binning strategy for calculating the cumulative
%                     histogram. 1 for using "histogram", 2 for using
%                     the data directly, and 3 for using "hist" with a
%                     specified number of bins.
%                     Optional. Default: 2.
%       numBinIn    : Number of bins to contruct histogram. 
%                     Required if binStrategy = 3, otherwise ignored.
%       plotName    : The title of the plotted figure.
%                     Optional. Default: 'Figure'.
%       meanLB      : Lower bound on mean of exponential.
%                     Optional. Default: 0.
%
%OUTPUT expParam    : Matrix with number of rows equal to number of fitted
%                     exponentials and 5 columns:
%                     Column 1: mean of each exponential
%                     Column 2: amplitude of each exponential
%                     Column 3: std of each mean
%                     Column 4: std of each amplitude
%                     Column 5: entry in the first row only storing the
%                     fit residuals variance.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS 
%
%Khuloud Jaqaman, June 2012

%% Output

expParam = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--fitHistWithExponentialsN: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%make sure observations is a column vector
observations = observations(:);

if nargin < 2 || isempty(alpha)
    alpha = 0.01;
else
    if alpha < 0 || alpha > 1
        disp('--fitHistWithExponentialsN: Variable "alpha" should be between 0 and 1!');
        errFlag = 1;
    end
end

if nargin < 3 || isempty(showPlot)
    showPlot = 1;
else
    if ~any(showPlot == [0,1,2])
        disp('--fitHistWithExponentialsN: Variable "showPlot" should be 0, 1, or 2!');
        errFlag = 1;
    end
end

if nargin < 4 || isempty(maxNumExp)
    maxNumExp = 100;
else
    if length(maxNumExp) > 1
        minNumExp = maxNumExp(1);
        maxNumExp = maxNumExp(2);
    else
        minNumExp = 1;
    end
    if minNumExp>1
        warning('minNumExp can only be taken into account if ''R''') %#ok<WNTAG>
    end
    if maxNumExp < 1
        disp('--fitHistWithExponentialsN: Variable "maxNumExp" should be at least 1!');
        errFlag = 1;
    end
end

if nargin < 5 || isempty(binStrategy)
    binStrategy = 2;
else
    if ~any(binStrategy == [1,2,3])
        disp('--fitHistWithExponentialsN: Variable "binStrategy" should be 1, 2 or 3!');
        errFlag = 1;
    end
end

if nargin < 6 || isempty(numBinIn)
    if binStrategy == 3
        disp('--fitHistWithExponentialsN: Please input value for numBinIn when binStrategy = 3!');
        errFlag = 1;        
    else
        numBinIn = [];
    end
end

if nargin < 7 || isempty(plotName)
    plotName = 'Figure';
end

if nargin < 8 || isempty(meanLB)
    meanLB = 0;
end

% check error flag
if errFlag
    return
end

%% Cumulative distribution function

%get the number of observations
numObservations = length(find(~isnan(observations)));

switch binStrategy
    
    case 1 %use "histogram"
        
        %calculate the histogram
        [numObsPerBin,binCenter] = optimalHistogram(observations,[],0);
        %         numObsPerBin = numObsPerBin'*(binCenter(2)-binCenter(1));
        numObsPerBin = numObsPerBin';
        binCenter = binCenter';
        
        %determine the number of bins used
        numBins = length(binCenter);
        
        %calculate the cumulative histogram
        cumHist = zeros(numBins,1);
        for iBin = 1 : numBins
            cumHist(iBin) = sum(numObsPerBin(1:iBin));
        end
        
    case 2
        
        %get the empirical cumulative distribution function
        [cumHist,binCenter] = ecdf(observations);
        binCenter = binCenter(2:end);
        cumHist = cumHist(2:end);
        numBins = length(binCenter);
        
        %downsample to about 1000 points if necessary
        if numBins > 1000
            dsIdx = unique(round(linspace(1,numBins,1000)))';
            cumHist = cumHist(dsIdx);
            binCenter = binCenter(dsIdx);
            numBins = length(dsIdx);
        end
        
        %make cumHist go from 1:numObservations
        cumHist = cumHist * numObservations;
        
    case 3
        
        %calculate the histogram
        [numObsPerBin,binCenter] = hist(observations,numBinIn);
        numObsPerBin = numObsPerBin';
        binCenter = binCenter';
        
        %determine the number of bins used
        numBins = length(binCenter);
        
        %calculate the cumulative histogram
        cumHist = zeros(numBins,1);
        for iBin = 1 : numBins
            cumHist(iBin) = sum(numObsPerBin(1:iBin));
        end
                
end

%% Fit

%initialize variables indicating number of fitted exponentials and their parameters
numExp = 0;
expParam = [];

%logical variable indicating whether to attempt to fit
fit = 1;

%set some optimization options
options = optimset('MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-10,'TolX',1e-10,'Display','off'); %,'Jacobian','on');

%maximum observed value, just to avoid repetition
maxVal = max(observations);

%fit the cumulative histogram with as many exponentials as necessary
while fit
    
    %add another exponential to the fit
    numExpT = numExp + 1;
    
    %calculate number of degrees of freedom
    numDegFreeT = numBins - 2*numExpT;
    
    %assign the parameter initial guesses
    meanIGInc = (0.9*maxVal-1.1*meanLB)/(numExpT+1);
    meanInitialGuess = (1.1*meanLB : meanIGInc : 0.9*maxVal)';
    meanInitialGuess = meanInitialGuess(2:end-1);
    ampInitialGuess = numObservations/numExpT*ones(numExpT,1);
    expParamT = [meanInitialGuess ampInitialGuess];
    
    %assign parameter initial values
    x0 = expParamT(:);
    
    %assign lower bounds
    lb = [meanLB*ones(numExpT,1) zeros(numExpT,1)];
    lb = lb(:);
    
    %assign upper bounds
    ub = [maxVal*ones(numExpT,1) 1.5*numObservations*ones(numExpT,1)];
    ub = ub(:);
    
    %estimate unknown parameters
    [param,~,residualsT,~,~,~,jacMatT] = lsqcurvefit(@calcCumDistrNExp,x0,...
        binCenter,cumHist,lb,ub,options);
    
    %get output from parameters vector
    expParamT = reshape(param,numExpT,2);
    
    %check whether addition of 1 exponential has significantly improved the fit
    if numExpT > 1 %if this is not the first fit
        
        %get test statistic, which is F-distributed
        testStat = (sum(residualsT.^2)/numDegFreeT)/...
            (sum(residuals.^2)/numDegFree);
        
        %get p-value of test statistic
        pValue = fcdf(testStat,numDegFree,numDegFreeT);
        
        %compare p-value to alpha
        %1-sided F-test: H0: F=1, H1: F<1
        if pValue <= alpha && numExpT <= maxNumExp %if p-value is smaller and the limit of exponentials isn't reached
            fit = 1; %accept this fit and attempt another fit with an additional exponential
        else %if p-value is larger
            fit = 0; %do not accept this fit and exit
        end
        
    end %(if numExpT > 1)
    
    %if this fit is accepted, update some variables
    if fit
        numExp = numExpT;
        expParam = expParamT;
        residuals = residualsT;
        numDegFree = numDegFreeT;
        jacMat = jacMatT;
    end
    
end %(while fit)

%calculate the parameter variance-covariance matrix
varCovMat = full((sum(residuals.^2)/numDegFree) * inv(jacMat'*jacMat));

%append standard deviation to output variable expParam
paramStd = sqrt(diag(varCovMat));
paramStd = reshape(paramStd,numExp,2);
expParam = [expParam paramStd];

%order the exponentials in ascending value of the mean
expMeans = expParam(:,1);
[~,orderIndx] = sort(expMeans);
expParam = expParam(orderIndx,:);

%get the p-values for the distances between adjacent exponentials and
%determine whether one exponential should be removed
if numExp > 1 && alpha < 1
    adjacentExpoDist = diff(expParam(:,1));
    distStd = sqrt( expParam(2:end,3).^2 + expParam(1:end-1,3).^2 );
    pValue = 1 - normcdf(adjacentExpoDist,0,distStd);
    if max(pValue) > 0.05
        removeExp = 1;
    else
        removeExp = 0;
    end
else
    removeExp = 0;
end

%remove exponentials if their means are not far enough from each other
while removeExp
    
    %remove one exponential from the fit
    numExp = numExp - 1;
    
    %calculate number of degrees of freedom
    numDegFree = numBins - 2*numExp;
    
    %assign the parameter initial guesses
    meanIGInc = (0.9*maxVal-1.1*meanLB)/(numExp+1);
    meanInitialGuess = (1.1*meanLB : meanIGInc : 0.9*maxVal)';
    meanInitialGuess = meanInitialGuess(2:end-1);
    ampInitialGuess = numObservations/numExp*ones(numExp,1);
    expParam = [meanInitialGuess ampInitialGuess];
    
    %assign parameter initial values
    x0 = expParam(:);
    
    %assign lower bounds
    lb = [meanLB*ones(numExp,1) zeros(numExp,1)];
    lb = lb(:);
    
    %assign upper bounds
    ub = [maxVal*ones(numExp,1) 1.5*numObservations*ones(numExp,1)];
    ub = ub(:);
    
    %estimate unknown parameters
    [param,~,residuals,~,~,~,jacMat] = lsqcurvefit(@calcCumDistrNExp,x0,...
        binCenter,cumHist,lb,ub,options);
    
    %get output from parameters vector
    expParam = reshape(param,numExp,2);
    
    %calculate the parameter variance-covariance matrix
    varCovMat = full((sum(residuals.^2)/numDegFree) * inv(jacMat'*jacMat));
    
    %append standard deviation to output variable expParam
    paramStd = sqrt(diag(varCovMat));
    paramStd = reshape(paramStd,numExp,2);
    expParam = [expParam paramStd];
    
    %order the exponentials in ascending value of the mean
    expMeans = expParam(:,1);
    [~,orderIndx] = sort(expMeans);
    expParam = expParam(orderIndx,:);
    
    %get the p-values for the distances between adjacent exponentials and
    %determine whether one exponential should be removed
    if numExp > 1
        adjacentExpoDist = diff(expParam(:,1));
        distStd = sqrt( expParam(2:end,3).^2 + expParam(1:end-1,3).^2 );
        pValue = 1 - normcdf(adjacentExpoDist,0,distStd);
        if max(pValue) > 0.05
            removeExp = 1;
        else
            removeExp = 0;
        end
    else
        removeExp = 0;
    end
    
end

%append the sum of squared residuals / # degrees of freedom to
%the first row of expParam
expParam(1,end+1) = sum(residuals.^2) / numDegFree;

%also append the BIC
expParam(1,end+1) = numBins*log(sum(residuals.^2)/numBins) + 2*numExpT*log(numBins);

%% Plotting

% make a histogram for plotting and output. Choose how to calculate bins
if showPlot == 2
    [numObsPerBinP,binCenterP] = optimalHistogram(observations,'smooth');
    numObsPerBinP = numObsPerBinP*(binCenterP(2)-binCenterP(1));
elseif showPlot ~= 0
    [numObsPerBinP,binCenterP] = optimalHistogram(observations);
    numObsPerBinP = numObsPerBinP*(binCenterP(2)-binCenterP(1));
end

%if the user wants to plot
if showPlot
    
    %get the distribution from the optimized parameters
    distrIndExp = zeros(numExp,length(binCenterP));
    for i=1:numExp
        distrIndExp(i,:) = expParam(i,2)*exppdf(binCenterP,expParam(i,1))*(binCenterP(2)-binCenterP(1));
    end
    distrNExp = sum(distrIndExp,1);
    
    %get the cumulative distribution from the optimized parameters
    cumDistrNExp = zeros(size(binCenter));
    for i=1:numExp
        cumDistrNExp = cumDistrNExp + expParam(i,2)*expcdf(binCenter,expParam(i,1));
    end
    
    %make new figure
    if isempty(plotName)
        figure
    else
        figure('Name',plotName,'NumberTitle','off')
    end
    
    %plot the histogram and the fitted exponentials in the left half of the
    %figure. Correct by the number of NaNs
    subplot(1,2,1);
    bar(binCenterP,numObsPerBinP,'k')
    hold on
    for i = 1 : numExp
        plot(binCenterP,distrIndExp(i,:) * sum(isfinite(observations))/...
            numObservations,'g--','LineWidth',2.5)
    end
    plot(binCenterP,distrNExp * sum(isfinite(observations))/...
        numObservations,'r','LineWidth',2.5)
    xlabel('Observation values')
    ylabel('Counts')
    
    %plot the cumulative histogram and the fitted exponentials in the right
    %half of the figure
    subplot(1,2,2);
    plot(binCenter,cumHist,'k.')
    hold on
    plot(binCenter,cumDistrNExp,'r','LineWidth',2.5)
    xlabel('Observation values')
    ylabel('Cumulative counts')

end %(if showPlot)

%%%%% ~~ the end ~~ %%%%%
