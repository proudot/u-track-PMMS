function testMaxCurvatureBezier()

clc;
clear all;

% Set some random control points
cP = [rand(4,2)*4  zeros(4,1)];
t0 = 0;
t1 = 1;

% Plot curve
curvePoints = renderBezier(cP,linspace(t0-0.5,t1+0.5,300)');
figure(1);
plot(curvePoints(:,1),curvePoints(:,2),'-r','LineWidth',3);
axis image;
axis ([-1 5 -1 5]);
axis off;
hold on;
curvePoints = renderBezier(cP,linspace(t0,t1,300)');
plot(curvePoints(:,1),curvePoints(:,2),'-b','LineWidth',1);

% Compute the maximum curvature
tic;
[curvature,t] = maxCurvatureBezier(cP,t0,t1);
toc;

% Print results
fprintf('The maximum curvature is %f\n',curvature);
fprintf('The maximum curvature is at t=%f\n',t);

% Plot the position of the maximum curvature
maxCurvaturePoint = renderBezier(cP,t);
scatter(maxCurvaturePoint(:,1),maxCurvaturePoint(:,2),'ob');

hold off;

end

