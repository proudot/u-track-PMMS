function distancePointBezierTest()

clc;
clear all;

% Set some random control points and draw the curve
cP = [[1 1 0];[4 2 0];[2 4 0];[4 4 0]];
% cP = [[0 0 0];[4 2 0];[2 4 0];[4 4 0];[3 0 0];[-3 6 0];[3 0 0]];

curvePoints = renderBezier(cP,linspace(0,1,300)');

figure(1);
plot(curvePoints(:,1),curvePoints(:,2),'-r','LineWidth',3);
axis image;
axis ([-1 5 -1 5]);
axis off;
hold on;

% Generate random points
nPoints = 4000; % 100
p = rand(nPoints,3)*5.5-0.75;
scatter(p(:,1),p(:,2),'.b');

% Compute the distance and draw the projections
for i=1:nPoints
    [d,t] = distancePointBezier(cP,p(i,:));
    projectionPoints = renderBezier(cP,t);
    plot([p(i,1),projectionPoints(1)],[p(i,2),projectionPoints(2)],'-k');
end

hold off;
end

