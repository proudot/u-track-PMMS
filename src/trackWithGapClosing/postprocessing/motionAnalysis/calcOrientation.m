function orientVec = calcOrientation(movieInfo,doPlot)
%CALCORIENTATION determines the direction of maximal scatter for a set of x and y coordinates
%
%SYNPOSIS orientVec = calcOrientation(movieInfo,doPlot)
%
%INPUT  movieInfo   : Output of detectSubResFeatures2D_StandAlone.
%       doPlot      : 1 to plot the x- and y-coordinates and the
%                     calculated direction of maximal scatter, 0 otherwise.
%                     Optional. Default: 1.
%
%OUTPUT orientVec   : Vector of direction of maximal scatter.
%
%Khuloud Jaqaman, December 2010

%% Output

orientVec = [];

%% Input

if nargin < 1
    disp('--calcOrientVec: Incorrect number of input arguments!');
    return
end
    
if nargin < 2 || isempty(doPlot)
    doPlot = 1;
end

%% Calculation

%get the x and y coordinates
xCoord = vertcat(movieInfo.xCoord);
yCoord = vertcat(movieInfo.yCoord);
xyCoord = [xCoord(:,1) yCoord(:,1)];

%calculate the coordinate covariance matrix
coordCovMat = cov(xyCoord);

%do eigenvalue decomposition to get the directions of minimal and maximal
%scatter
[eigVec,eigVal] = eig(coordCovMat);
eigVal = diag(eigVal);

%get the direction of maximal scatter
orientVec = eigVec(:,eigVal==max(eigVal));

%% Plot

if doPlot
    
    %calculate center for plotting
    centerCoord = mean(xyCoord);
    
    %plot
    figure, hold on
    plot(xyCoord(:,1),xyCoord(:,2),'.')
    plot([centerCoord(:,1) centerCoord(:,1)+100*orientVec(1)],...
        [centerCoord(:,2) centerCoord(:,2)+100*orientVec(2)],'r')
    
end

