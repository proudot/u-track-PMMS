function [T,normalT] = tangentBezier(P, t)
% function [T,normalT] = tangentBezier(P, t)
% tangentBezier computes the unnormalized and normalized tangent vector at 
% the positions defined by the node vector t.
%
% Required Inputs:
% P              A n+1 x d array representing the d-dimensional control 
%                points of the Bezier curve of degree n.
%
% t              1 x m array representing the nodes at which the tangent
%                should be evaluated.
% 
% Optional Inputs:
%
% Outputs:
% T              A m x d array representing the d-dimensional
%                unnormalized tangent vectors at the nodes defined by t
%
% normalT        A m x d array representing the d-dimensional
%                normalized tangent vectors at the nodes defined by t
%
% Pascal Berard, October 2011

% Determine the curve degree and the number of nodes
n = size(P,1)-1;
m = numel(t);

% Compute the Benrstein Matrix of order n-1
Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
Bn_1 = (bsxfun(@power, t, 0:n-1) .* bsxfun(@power, 1 - t, n-1:-1:0)) * Cn_1k;

% Compute the derivative of the Bernstein matrix against t
z = zeros(m,1);
Bt = n * ([z Bn_1] - [Bn_1 z]);

% Compute the tangent vector
T = Bt * P;

% Compute the normalized tangent vector if needed
if nargout > 1
    dim = size(P,2);
    normalT = T./repmat(sqrt(sum(T.^2,2)),1,dim);
end
