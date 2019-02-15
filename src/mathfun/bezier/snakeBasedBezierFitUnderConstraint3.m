function [C,T,res,lambda] = snakeBasedBezierFitUnderConstraint3(data,n,beta,varargin)
% snakeBasedBezierFit computes the optimal n^th Bezier curves that fits
% the data points using a scale-invariant snake-based functional.
%
% Required Inputs:
% data           A m x d array representing a set of d-dimensional points
% n              Degree of the Bezier curve.
% beta           Regularization weight parameter.
%
% Optional Inputs:
% SigmaX         A d x d matrix which represents the variance-covariance
%                matrix of the localization precision of the data.
% MaxFunEvals    Maximum number of fonctional evaluations during lsqnonlin.
% MaxIter        Maximum number of interations during lsqnonlin.
% Display        Verbose mode during lsqnonlin.
% TolX           Tolerance on the solution (i.e. t) during lsqnonlin.
% TolFun         Tolerance on the functional during lsqnonlin.
%
% Outputs:
%
% C:         Control points of the optimal Bezier curve. a (n+1)xd matrix
% T:         a mx1 vector
% res:       a mx3 residual vector corresponding to the distances
%            between data points and the curve. This is not taking into account
%            SigmaX.
% lambda     a structure lambda whose fields contain the Lagrange
%            multipliers at the solution [C;T]
%
% Sylvain Berlemont, Dec. 2011

%% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(data) isnumeric(data) && size(data,2) > 1);
ip.addRequired('n', @(n) n > 0);
ip.addRequired('beta', @(beta) beta >= 0);
ip.addParamValue('SigmaX', eye(size(data,2)), @isnumeric);
ip.addParamValue('Algorithm', 'interior-point', @isstr);
ip.addParamValue('MaxFunEvals', 1e4, @isscalar);
ip.addParamValue('MaxIter', 1e4, @isscalar);
ip.addParamValue('Display', 'off', @isstr);
ip.addParamValue('TolX', 1e-8, @isscalar);
ip.addParamValue('TolFun', 1e-8, @isscalar);

ip.parse(data, n, beta, varargin{:});
SigmaX = ip.Results.SigmaX;
algorithm = ip.Results.Algorithm;
maxFunEvals = ip.Results.MaxFunEvals;
maxIter = ip.Results.MaxIter;
display = ip.Results.Display;
tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

%% Setup the optimization algorithm
opts = optimset('Algorithm', algorithm, ...
  'DerivativeCheck', 'off', ...
  'Display', display, ...
  'FunValCheck', 'off', ...
  'GradObj', 'off', ...
  'GradConstr', 'on', ...
  'Hessian', 'off', ...
  'MaxFunEvals', maxFunEvals, ...
  'MaxIter', maxIter, ...
  'TolX', tolX, ...
  'Tolfun', tolFun);

% Array of function handlers for regularization term computation
regFuncs = {@computeRegTermN1D2, @computeRegTermN2D2, @computeRegTermN3D2;...
  @computeRegTermN1D3, @computeRegTermN2D3, @computeRegTermN3D3};

%% Compute an initial solution
[m d] = size(data);

% Weight beta
beta = beta * m;

% Number of parameters of the fitting problem is equal to
% - number of control point coordinates: d * (n+1)
% - number of nodes without the first and last ones: m-2
% pDim = d * (n+1) + (m-2);

% Compute the initial nodes
T = linspace(0,1,m)';
% Compute the initial control points
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

B = (bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0)) * Cnk;
[Q1 R11] = qr(B,0);
C = R11 \ (Q1' * data);
% Note we do not optimize t1 and tm since t1=0 and tm=1 in all cases.
X = [C(:);T(2:end-1)];

% Constraints for T value in [0,1]
lb = -inf(size(X));
ub = +inf(size(X));
lb(d * (n+1) + 1:end) = 0;
ub(d * (n+1) + 1:end) = 1;

%% Optimization
[X,~,~,~,lambda] = fmincon(@fun,X,[],[],[],[],lb,ub,[],opts);

% Compute the residual
T = [0; X(d * (n + 1) + 1:end); 1];
B = (bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0)) * Cnk;
res = B * C - data;

  function F = fun(X)
    
    % Retrieve the control point coordinates from X
    C = reshape(X(1:d * (n + 1)), n + 1, d);
    
    % Retrieve the nodes from X and add t0 and tm
    T = [0; X(d * (n + 1) + 1:end); 1];
    
    % Compute the Bernstein matrix
    B = (bsxfun(@power, T, 0:n) .* bsxfun(@power, 1 - T, n:-1:0)) * Cnk;
    
    % Compute the data fidelity term
    V = data - B * C;
    res = sum((V / SigmaX) .* V, 2);
    dataFidelity = sum(res);
    
    % Append the regularization terms and the contraints
    F = dataFidelity + beta * regFuncs{d-1,n}(C);
  end
end

function reg = computeRegTermN1D2(~)

reg = 0;

end

function reg = computeRegTermN1D3(~)

reg = 0;

end

function reg = computeRegTermN2D2(C)

CC = num2cell(C);
[x0, x1, x2, y0, y1, y2] = CC{:};

ds2 = 4 * ((x0 - 2 * x1 + x2)^2 + (y0 - 2 * y1 + y2)^2);

l2 = (4/3) * (x0^2 + x1^2 - x1 * x2 + x2^2 - x0 * (x1 + x2) + ...
    y0^2 - y0 * y1 + y1^2 - (y0 + y1) * y2 + y2^2);

reg = ds2 / l2;

end

function reg = computeRegTermN2D3(C)

CC = num2cell(C);
[x0, x1, x2, y0, y1, y2, z0, z1, z2] = CC{:};

ds2 = 4 * ((x0 - 2 * x1 + x2)^2 + (y0 - 2 * y1 + y2)^2 + (z0 - 2 * z1 + z2)^2);

l2 = (4/3) * (x0^2 + x1^2 - x1 * x2 + x2^2 - x0 * (x1 + x2) + y0^2 - ...
    y0 * y1 + y1^2 - (y0 + y1) * y2 + y2^2 + z0^2 - z0 * z1 + z1^2 - ...
    (z0 + z1) * z2 + z2^2);

reg = ds2 / l2;

end

function reg = computeRegTermN3D2(C)

CC = num2cell(C);
[x0, x1, x2, x3, y0, y1, y2, y3] = CC{:};

ds2 = 12 * (x0^2 + 3 * (x1^2 - x1 * x2 + x2^2) - 3 * x2 * x3 + x3^2 + ... 
   x0 * (x3 - 3 * x1) + y0^2 + 3 * (y1^2 - y1 * y2 + y2^2) - 3 * y2 * y3 + ...
   y3^2 + y0 * (y3 - 3 * y1));

l2 = (3/5) * (3 * x0^2 + 2 * x1^2 + 2 * x2^2 + x1 * (x2 - 2 * x3) - ...
    3 * x2 * x3 + 3 * x3^2 - x0 * (3 * x1 + 2 * x2 + x3) + 3 * y0^2 - ...
    3 * y0 * y1 + 2 * y1^2 - 2 * y0 * y2 + y1 * y2 + 2 * y2^2 - ...
    (y0 + 2 * y1 + 3 * y2) * y3 + 3 * y3^2);

reg = ds2 / l2;

end

function reg = computeRegTermN3D3(C)

CC = num2cell(C);
[x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3] = CC{:};

ds2 = 12 * (x0^2 + 3 * (x1^2 - x1 * x2 + x2^2) - 3 * x2 * x3 + x3^2 + ...
   x0 * (x3 -3 * x1) + y0^2 + 3 * (y1^2 - y1 * y2 + y2^2) - ...
   3 * y2 * y3 + y3^2 + y0 * (y3 -3 * y1) + z0^2 + 3 * (z1^2 - ...
   z1 * z2 + z2^2) - 3 * z2 * z3 + z3^2 + z0 * (z3 -3 * z1));

l2 = (3/5) * (3 * x0^2 + 2 * x1^2 + 2 * x2^2 + x1 * (x2 - 2 * x3) - ...
    3 * x2 * x3 + 3 * x3^2 - x0 * (3 * x1 + 2 * x2 + x3) + 3 * y0^2 - ...
    3 * y0 * y1 + 2 * y1^2 - 2 * y0 * y2 + y1 * y2 + 2 * y2^2 - ...
    (y0 + 2 * y1 + 3 * y2) * y3 + 3 * y3^2 + 3 * z0^2 - 3 * z0 * z1 + ...
    2 * z1^2 - 2 * z0 * z2 + z1 * z2 + 2 * z2^2 - (z0 + 2 * z1 + 3 * z2) * ...
    z3 + 3 * z3^2);

reg = ds2 / l2;

end