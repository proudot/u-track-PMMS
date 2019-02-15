function [cP,t] = truncateBezier(cP,tStart,tEnd,t)
% function [cP,t] = truncateBezier(cP,tStart,tEnd,t)
% truncateBezier computes the control points of a interval of the input 
% Bezier curve. Optionally, points of the old curve can be mapped onto the
% new curve.
%
% Required Inputs:
% cP                A N x D array representing a set of 3-dimensional 
%                   control points where N is 2,3 or 4. D is the dimension
%                   of the control points.
% tStart            Defines the start point of the segment            
% tEnd              Defines the end point of the segment
% 
% Optional Inputs:
% t                 M x 1 array containing parameter values of the old
%                   curve that will be mapped to the new one
%
% Outputs:
% cP                The control points of the new curve
% t                 The transformed parametrization 
%
% Pascal Berard, March 2012

% Transform the parametrization
if nargin == 4
    t = (t-tStart)/(tEnd-tStart);
end

% Determine the curve degree
n = size(cP, 1) - 1;

% Sample curve parameter
tSample = linspace(0,1,n+1)';
samplePoints = renderBezier(cP,tSample);

% Transform the parametrization
tSample = (tSample-tStart)/(tEnd-tStart);

% Compute the Bernstein Matrix
B = bsxfun(@power, tSample, 0:n) .* bsxfun(@power, 1 - tSample, n:-1:0);
B = B * diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);

% Compute the new control points
cP = B\samplePoints;

