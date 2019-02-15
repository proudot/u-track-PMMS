function testComputeICP
% This function test the Iterative Closest Point algorithm (ICP).
%
% See computeICP.m for details.
%
% Sylvain Berlemont, Dec. 2009

tol = 1e-9; % ending condition of ICP
prec = 10; % number of significant digit
maxIter = 10;

% Test 1.
X1 = rand(100, 3);

[T R] = computeICP(X1, X1, maxIter, tol);
T = round(T * 10^prec) / 10^prec;
R = round(R * 10^prec) / 10^prec;
if all(T == 0) && all(all(R == eye(3)))
    disp('Test 1: OK');
else
    disp('Test 1: NOT OK');
end

% Test 2
X1 = rand(100, 3);
X2 = X1;
shift = [1e-9 0 0];
X1 = bsxfun(@plus, X1, shift);

[T R] = computeICP(X1, X2, maxIter, tol);
T = round(T * 10^prec) / 10^prec;
R = round(R * 10^prec) / 10^prec;
if all(T == shift') && all(all(R == eye(3)))
    disp('Test 2: OK');
else
    disp('Test 2: NOT OK');
end

% Test 3
X1 = rand(100, 3);
theta = (.5 - rand(1)) * pi;
ct = cos(theta);
st = sin(theta);
Rref = [ct st 0; -st ct 0; 0 0 1];
X2 = X1;
c1 = mean(X1);
X1 = bsxfun(@minus, X1, c1);
for i = 1:size(X1, 1)
    X1(i, :) = X1(i, :) * Rref;
end
X1 = bsxfun(@plus, X1, c1);

[T R] = computeICP(X1, X2, maxIter, tol);
T = round(T * 10^prec) / 10^prec;
R = round(R * 10^prec) / 10^prec;
if all(T == 0) && abs(acos(R(1))^2 - theta^2) < 1e-9 && R(end) == 1
    disp('Test 3: OK');
else
    disp('Test 3: NOT OK');
end

% Test 4
X1 = rand(100, 3);
shift = [1e-9 0 0];
theta = (.5 - rand(1)) * pi;
ct = cos(theta);
st = sin(theta);
Rref = [ct st 0; -st ct 0; 0 0 1];
X2 = X1;
c1 = mean(X1);
X1 = bsxfun(@minus, X1, c1);
for i = 1:size(X1, 1)
    X1(i, :) = X1(i, :) * Rref;
end
X1 = bsxfun(@plus, X1, c1 + shift);

[T R] = computeICP(X1, X2, maxIter, tol);
T = round(T * 10^prec) / 10^prec;
R = round(R * 10^prec) / 10^prec;
if all(T == shift') && abs(acos(R(1))^2 - theta^2) < 1e-9 && R(end) == 1
    disp('Test 4: OK');
else
    disp('Test 4: NOT OK');
end
