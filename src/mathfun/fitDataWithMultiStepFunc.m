function [stepX,valY] = fitDataWithMultiStepFunc(x,y,numStepsRange,alpha,doPlot)
%FITDATAWITHMULTISTEPFUNC fits a time series with a multi-step function with number of steps determined automatically
%
%SYNOPSIS [stepX,valY] = fitDataWithMultiStepFunc(x,y,numStepsRange,alpha,doPlot)
%
%INPUT  x            : Independent variable of time series.
%       y            : Dependent variable of time series.
%       numStepsRange: Row vector indicating minimum and maximum number
%                      steps to attempt to fit.
%                      Optional. Default: [0 5].
%       alpha        : Alpha-value for F-test to determine number of steps.
%                      Optional. Default: 0.05.
%       doPlot       : 1 to plot data and fit, 0 otherwise.
%                      Optional. Default: 0.
%
%OUTPUT stepX        : x-values at which there is a step.
%       valY         : y-values between steps. If there are n steps, there
%                      will be n+1 y-values. First value is before first
%                      step, etc., and last value is after last step.
%
%Khuloud Jaqaman, August 2014

%% Output
stepX = [];
valY = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--fitDataWithMultiStepFunc: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(numStepsRange)
    numStepsRange = [0 5];
end

if nargin < 4 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end

%get range of steps to try to fit
minNumSteps = numStepsRange(1);
maxNumSteps = numStepsRange(2);

%determine number of initial guesses for fitting
numGuess = max(10,maxNumSteps*2);

%take moving average of data to help with step initialization
[~,yMA] = movAvSmooth(y,5,2,0);

%get rid of NaNs in data
indx = find(~isnan(y));
x = x(indx);
y = y(indx);
yMA = yMA(indx);
numDataPoints = length(indx);

%% Calculation

%locate positions of largest jumps in data for step initialization
yMAdiff = abs(diff(yMA));
[~,iSorted] = sort(yMAdiff,1,'descend');
stepPosGuess = x(iSorted);
numGuess = min(numGuess,length(stepPosGuess));

%minimization options
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-3);

%determine parameters of "0 steps" fit, which is nothing but data mean
%this is the minimal fit
param = mean(y);
residuals = param - y;
numDegFree = numDataPoints - 1;

%fit larger models one by one, checking residuals
fit = 1;
iStep = 0;
while fit && iStep < maxNumSteps
    
    %number of steps for current fit
    iStep = iStep + 1;
    
    %number of parameters and degrees of freedom for current fit
    numParamT = 2*iStep + 1;
    numDegFreeT = numDataPoints - numParamT;
    
    %intermediate variables to handle multiple initial guesses
    paramT = zeros(numParamT,numGuess);
    residualsT = zeros(numDataPoints,numGuess);
    ssr = zeros(numGuess,1);
    
    %use different initial guesses to increase chance of finding global
    %minimum
    for iGuess = 1 : numGuess
        
        %parameter initial guess
        stepX0 = sort([param(1:iStep-1); stepPosGuess(iGuess)]);
        stepIndx = zeros(iStep,1);
        for iTmp = 1 : iStep
            stepIndx(iTmp) = find(abs(x-stepX0(iTmp))==min(abs(x-stepX0(iTmp))));
        end
        stepIndx = [1; stepIndx; numDataPoints]; %#ok<AGROW>
        valY0 = ones(iStep+1,1);
        for iPart = 1 : iStep+1
            valY0(iPart) = mean(y(stepIndx(iPart):stepIndx(iPart+1)));
        end
        
        %fit
        paramT(:,iGuess) = fminsearch(@multiStepFunction,[stepX0; valY0],options,x,y);
        [~,residualsT(:,iGuess)] = multiStepFunction(paramT(:,iGuess),x,y);
        ssr(iGuess) = sum(residualsT(:,iGuess).^2);
        
    end
    
    %keep fit with smallest SSR
    indxGood = find(ssr == min(ssr));
    indxGood = indxGood(1);
    paramT = paramT(:,indxGood);
    residualsT = residualsT(:,indxGood);
    
    %keep fit if number of steps <= minimum, otherwise compare residuals
    %and retain only if fit is significantly better than previous fit
    if iStep <= minNumSteps
        param = paramT;
        numDegFree = numDegFreeT;
        residuals = residualsT;
    else
        testStat = (sum(residualsT.^2)/numDegFreeT)/...
            (sum(residuals.^2)/numDegFree);
        pValue = fcdf(testStat,numDegFree,numDegFreeT);
        if pValue < alpha
            param = paramT;
            numDegFree = numDegFreeT;
            residuals = residualsT;
        else
            fit = 0;
        end
    end
    
end

%parameters for output
numSteps = (length(param)-1)/2;
stepX = param(1:numSteps);
valY = param(numSteps+1:end);

%% Plotting

if doPlot
    
    overlayMultiStepFunc(x,y,stepX,valY)
    
%     figure, hold on
%     plot(x,y)
%     if numSteps == 0
%         plot([x(1) x(end)],valY*[1 1],'r');
%     else
%         plot([x(1) stepX(1)],valY(1)*[1 1],'r');
%         plot(stepX(1)*[1 1]',valY(1:2),'r');
%         for iStep = 1 : numSteps-1
%             plot([stepX(iStep) stepX(iStep+1)],valY(iStep+1)*[1 1],'r');
%             plot(stepX(iStep+1)*[1 1]',valY(iStep+1:iStep+2),'r');
%         end
%         plot([stepX(end) x(end)],valY(end)*[1 1],'r');
%     end
    
end


%% ~~~ the end ~~~

