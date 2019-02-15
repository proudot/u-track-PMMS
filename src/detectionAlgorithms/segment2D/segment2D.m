function F = segment2D(x, y, amp, l, sigma, theta, xRange, yRange, nzIdx)
% sub-resolution 2D segment model defined by 6 parameters:
%    xy      : position of the segment's center
%    amp     : mean amplitude along the segment
%    l       : length
%    sigma   : half width of the segment
%    theta   : orientation [-pi/2, pi/2)
%
% F = segment2D(x, y, amp, l, sigma, theta, xRange, yRange, nzIdx)
%
% parameters:
%
% (x,y)              position of the segment's center
%
% amp                amplitude of the segment
%
% l                  length of the segment
%
% sigma              half width of the segment
%
% theta              orientation of the segment
%
% (xRange, yRange)   2 vectors representing the 2-dimensional support of
%                    the segment. This support can be determined using
%                    segment2DSupport() function.
%
% nzIdx              linear indices of a NxM matrix (N = numel(yRange) and
%                    M = numel(xRange)) where the model is defined. If not
%                    provided, nzIdx = 1:N*M. These indices can be
%                    determined using segment2DSupport() function.
%
% output:
% F                  the model defined on nzIdx pixels.
%
% Sylvain Berlemont, 2010

xRange = xRange - x;
yRange = yRange - y;

N = numel(yRange);
M = numel(xRange);

if nargin < 9 || isempty(nzIdx)
    nzIdx = 1:N*M;
end

[X Y] = meshgrid(xRange, yRange);
X = X(nzIdx);
Y = Y(nzIdx);

ct = cos(theta);
st = sin(theta);
    
C0 = .5*2^(-1/2)*sigma^(-1);
C1 = exp(-.5 * sigma^(-2) * (Y * ct - X * st).^2);
C2 = erf(C0 * (l + 2 * X * ct + 2 * Y * st));
C3 = erf(C0 * (l - 2 * X * ct - 2 * Y * st));
C6 = erf(l * C0);

F = .5 * amp * C1 .* (C2 + C3) ./ C6;
