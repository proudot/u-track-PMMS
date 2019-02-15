function overlayMultiStepFunc(x,y,stepX,valY,figTitle)
%OVERLAYMULTISTEPFUNC overlays a multi-step function on time series
%
%SYNOPSIS overlayMultiStepFunc(x,y,stepX,valY)
%
%INPUT  x            : Independent variable of time series.
%       y            : Dependent variable of time series.
%       stepX        : x-values at which there is a step.
%       valY         : y-values between steps. If there are n steps, there
%                      will be n+1 y-values. First value is before first
%                      step, etc., and last value is after last step.
%       figTitle     : Figure title.
%
%OUTPUT the plot
%
%Khuloud Jaqaman, August 2014

if nargin < 5 || isempty(figTitle)
    figure
else
    figure('Name',figTitle,'NumberTitle','off')
end
hold on

plot(x,y,'marker','.')

numSteps = length(stepX);
if numSteps == 0
    plot([x(1) x(end)],valY*[1 1],'r');
else
    plot([x(1) stepX(1)],valY(1)*[1 1],'r');
    plot(stepX(1)*[1 1]',valY(1:2),'r');
    for iStep = 1 : numSteps-1
        plot([stepX(iStep) stepX(iStep+1)],valY(iStep+1)*[1 1],'r');
        plot(stepX(iStep+1)*[1 1]',valY(iStep+1:iStep+2),'r');
    end
    plot([stepX(end) x(end)],valY(end)*[1 1],'r');
end


%% ~~~ the end ~~~

