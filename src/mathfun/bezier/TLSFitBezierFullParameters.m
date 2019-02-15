function [P t res] = TLSFitBezierFullParameters(data, n, varargin)
% TLSFitBezier computes the best total least squares fit of a nth-degree
% Bezier curve to a set of ordered data points. As oppose to the function
% TLSFitBezier, the optimization is performed on the variable t1, ..., tm,
% x0, ..., xn, y0, ... yn.
%
% Required Inputs:
% data           A m x d array representing a set of d-dimensional points
% n              Degree of the Bezier curve.
% 
% Optional Inputs:
% MaxFunEvals    Maximum number of fonctional evaluations during lsqnonlin.
% MaxIter        Maximum number of interations during lsqnonlin.
% Display        Verbose mode during lsqnonlin.
% TolX           Tolerance on the solution (i.e. t) during lsqnonlin.
% TolFun         Tolerance on the functional during lsqnonlin.
%
% Outputs:
% P              A n+1 x d array representing the set of d-dimensional
%                control points defining the Bezier curve.
%
% t              A m x 1 array of parameter value for each input points. t
%                belongs to [0, 1].
%
% res            A m x 1 array of residue, i.e the orthogonal distance
%                between a data point and the fitted curve.
%
% Sylvain Berlemont, August 2011

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isnumeric);
ip.addRequired('n', @(n) n > 0);
ip.addParamValue('MaxFunEvals', 1e4, @isscalar);
ip.addParamValue('MaxIter', 1e4, @isscalar);
ip.addParamValue('Display', 'on', @isstring);
ip.addParamValue('TolX', 1e-8, @isscalar);
ip.addParamValue('TolFun', 1e-8, @isscalar);

ip.parse(data, n, varargin{:});
maxFunEvals = ip.Results.MaxFunEvals;
maxIter = ip.Results.MaxIter;
display = ip.Results.Display;
tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

% Define options of lsqnonlin algorithm
opts = optimset('Jacobian', 'off', ...
  'MaxFunEvals', maxFunEvals, ...
  'MaxIter', maxIter, ...
  'Display', display, ...
  'TolX', tolX, ...
  'Tolfun', tolFun);

% Define initial vector of nodes t0
[m dim] = size(data);
t = linspace(0,1,m)';

% Compute an initial set of control points
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;
[Q1 R11] = qr(B,0);
P = R11 \ (Q1' * data);

% Define the initial set of parameters
X = [t; P(:)];

% Solve the non-linear optimization on t
Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
fun = @(X) r(X, data, Cnk, Cn_1k, n);
lb = [zeros(size(t)); -inf(size(P(:)))];
ub = [ones(size(t)); +inf(size(P(:)))];
[X, ~, res] = lsqnonlin(fun, X, lb, ub, opts);

t = X(1:m);

P = reshape(X(m+1:end), [n+1, dim]);

% Reshape residual
res = sqrt(sum(reshape(res, [m, dim]).^2, 2));

function [F J] = r(X, data, Cnk, Cn_1k, n)
    
[m dim] = size(data);

t = X(1:m);
P = reshape(X(m+1:end), [n+1,dim]);

% Compute Bernstein Matrix
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

% Compute the residual
F = data - B * P;
F = F(:);

if nargout > 1
  % Compute the Benstein Matrix of order n-1
  Bn_1 = (bsxfun(@power, t, 0:n-1) .* bsxfun(@power, 1 - t, n-1:-1:0)) * Cn_1k;
  
  % Compute derivative of B against t
  z = zeros(m,1);
  Bt = n * ([z Bn_1] - [Bn_1 z]);
  
  J11 = diag(-Bt * P(:,1));
  J21 = diag(-Bt * P(:,2));
  
  z = zeros(m, n+1);
  J = [J11, -B, z; J21, z, -B];
end