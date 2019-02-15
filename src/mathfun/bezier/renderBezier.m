function Y = renderBezier(P, t)
% function Y = renderBezier(P, t)
% renderBezier evaluates the Bezier curve at the positions defined by the 
% node vector t.
%
% Required Inputs:
% P              A n+1 x d array representing the d-dimensional control 
%                points of the Bezier curve of degree n.
%
% t              1 x m or m x 1 array representing the nodes at which the curve
%                should be evaluated.
% 
% Optional Inputs:
%
% Outputs:
% Y              A m x d array representing the d-dimensional
%                unnormalized tangent vectors at the nodes defined by t
%
% Pascal Berard, October 2011

% Reshape t
t = reshape(t,numel(t),1);

% Determine the curve degree
n = size(P, 1) - 1;

% Compute the Bernstein Matrix
B = bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0);
B = B * diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

% Compute the curve points
Y = B * P;
