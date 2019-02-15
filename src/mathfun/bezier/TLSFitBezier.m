function [P t res] = TLSFitBezier(X, n, varargin)
% TLSFitBezier computes the best total least squares fit of a nth-degree
% Bezier curve to a set of ordered data points.
%
% Required Inputs:
% X              A m x d array representing a set of d-dimensional points
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
% Sylvain Berlemont, May 2011

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('X', @isnumeric);
ip.addRequired('n', @(n) n > 0);
ip.addParamValue('MaxFunEvals', 1e4, @isscalar);
ip.addParamValue('MaxIter', 1e4, @isscalar);
ip.addParamValue('Display', 'off', @isstring);
ip.addParamValue('TolX', 1e-8, @isscalar);
ip.addParamValue('TolFun', 1e-8, @isscalar);

ip.parse(X, n, varargin{:});
maxFunEvals = ip.Results.MaxFunEvals;
maxIter = ip.Results.MaxIter;
display = ip.Results.Display;
tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

% Define options of lsqnonlin algorithm
opts = optimset('Jacobian', 'on', ...
  'MaxFunEvals', maxFunEvals, ...
  'MaxIter', maxIter, ...
  'Display', display, ...
  'TolX', tolX, ...
  'Tolfun', tolFun);

% Define initial vector of nodes t0
[m d] = size(X);
t = linspace(0,1,m)';

% Solve the non-linear optimization on t
Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
fun = @(t) r(t, X, Cnk, Cn_1k, n);
[t, ~, res] = lsqnonlin(fun, t, zeros(size(t)), ones(size(t)), opts);

% Compute the control points
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;
[Q1 R11] = qr(B,0);
P = R11 \ (Q1' * X);

% Reshape residual
res = sqrt(sum(reshape(res, [m, d]).^2, 2));

function [F J] = r(t, X, Cnk, Cn_1k, n)
    
[m dim] = size(X);

% Compute Bernstein Matrix
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

% Compute residual = (Id - Bn * pinv(Bn)) * X
[Q1 R11 EE] = qr(B,0);
Q2Q2t = eye(m) - Q1 * Q1';
r = Q2Q2t * X;
F = r(:);

if nargout > 1
  % Compute the Benstein Matrix of order n-1
  Bn_1 = (bsxfun(@power, t, 0:n-1) .* bsxfun(@power, 1 - t, n-1:-1:0)) * Cn_1k;
  
  % Compute derivative of B against t
  z = zeros(m,1);
  Bt = n * ([z Bn_1] - [Bn_1 z]);
  
  % Compute P
  E = zeros(n + 1);
  E(sub2ind(size(E), EE, 1:(n+1))) = 1;
  P = Bt * E * (R11 \ Q1');

  % Compute Jacobian Matrix
  J = zeros(m * dim, m);
  for d = 1:dim
    J((d-1) * m + 1:d * m,:) = -(Q2Q2t * diag(P * X(:,d)) + P' * diag(r(:,d)));
  end
end