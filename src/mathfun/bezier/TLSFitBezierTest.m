function TLSFitBezierTest(n)

m = 60;
x = linspace(0,2 * pi, m);
y = 1*sin(x) + rand(size(x)) * .5;
data = [x' y'];

figure(11);
plot(x,y,'ro');

tic;
% [P1, t1] = TLSFitBezier([x' y'], n);
% [P2, t2] = TLSFitBezierFullParameters([x' y'], n);
[P3, t3] = TLSFitBezierConstraint1([x' y'], n);
% w = ones(size(data));
% w = [1000*ones(size(data,1),1),1000*ones(size(data,1),1)];
% [P3, t3] = TLSFitBezierWeightedConstrainedCP([x' y'], w, n); % Does not yet work in 2D!
toc;

% C1 = renderBezier(P1, t1);
% C2 = renderBezier(P2, t2);
C3 = renderBezier(P3, sort(t3));
C3_unordered = renderBezier(P3, t3);
% C3 = renderBezier(P3, linspace(0,1,200)');

hold on;

% Plot nodes and control points
% plot(C1(:,1),C1(:,2),'b-');
% plot(C2(:,1),C2(:,2),'r-');
plot(C3(:,1),C3(:,2),'g.-');
plot(P3(:,1),P3(:,2),'bo');

% Plot residuals
for p=1:m
    plot([x(p) C3_unordered(p,1)],[y(p) C3_unordered(p,2)],'b-');
end

axis equal;
hold off;

% % Compute the dot-product of the first data point and control point
% u = [x(1) - P1(1,1), y(1) - P1(1,2)];
% v = [P1(2, 1) - P1(1, 1), P1(2, 2) - P1(1, 2)];
% dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
% fprintf(1, 'dot product at start point (unconstraint) = %E\n', dot);
% 
% u = [x(m) - P1(n+1,1), y(m) - P1(n+1,2)];
% v = [P1(n, 1) - P1(n+1, 1), P1(n, 2) - P1(n+1, 2)];
% dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
% fprintf(1, 'dot product at end point (unconstraint) = %E\n', dot);
% 
% u = [x(1) - P2(1,1), y(1) - P2(1,2)];
% v = [P2(2, 1) - P2(1, 1), P2(2, 2) - P2(1, 2)];
% dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
% fprintf(1, 'dot product at start point (full param) = %E\n', dot);
% 
% u = [x(m) - P2(n+1,1), y(m) - P2(n+1,2)];
% v = [P2(n+1, 1) - P2(n, 1), P1(n+1, 2) - P2(n, 2)];
% dot = sum(u .* v) ./ (sqrt(sum(u.^2)) * sqrt(sum(v.^2)));
% fprintf(1, 'dot product at end point (full param) = %E\n', dot);
