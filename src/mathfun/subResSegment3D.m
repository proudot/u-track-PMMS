function F = subResSegment3D(xRange, yRange, zRange, A, Bg, sigmaPSF, l, theta1, theta2, nzIdx)
% Sub-resolution Segment Model
% F = subResSegment3D(xRange, yRange, zRange, A, Bg, sigmaPSF, l, theta1, theta2, nzIdx)
%
% parameters:
% (xRange, yRange,zRange)   3 vectors representing the 3-dimensional
%                           support of the segment. This support can be
%                           determined using subResSegment3DSupport()
%                           function.
%
% A                         amplitude of the segment
%
% Bg                        amplitude of the background
%
% sigmaPSF                  half width of the gaussian PSF model.
%
% l                         length of the segment
%
% theta1                    orientation of the segment on xOy plan
%
% theta2                    orientation of the segment on x0z plan
%
% nzIdx                     linear indices of a NxMxL matrix (N =
%                           numel(yRange), M = numel(xRange) and L =
%                           numel(zRange)) where the model is defined. If
%                           not provided, nzIdx = 1:N*M*L. These indices
%                           can be determined using
%                           subResSegment3DSupport() function.
%
% output:
% F is a NxMxL matrix where N = numel(Y), M = numel(X) and L = numel(Z)
%
% Sylvain Berlemont, 2010

N = numel(yRange);
M = numel(xRange);
L = numel(zRange);

if nargin < 7 || isempty(nzIdx)
    nzIdx = 1:N*M*L;
end

[X Y Z] = meshgrid(xRange, yRange, zRange);
X = X(nzIdx);
Y = Y(nzIdx);
Z = Z(nzIdx);

ct1 = cos(theta1);
ct2 = cos(theta2);
st1 = sin(theta1);
st2 = sin(theta2);

l = l / 2;
c0 = sqrt(2) * sigmaPSF;
c = A / (2 * erf(l / c0));

tmp1 = l + Z * ct1 + X * ct2 * st1 + Y * st1 * st2;
tmp2 = l - Z * ct1 + X * ct2 * st1 - Y * st1 * st2;

F = zeros(N,M,L);

F(nzIdx) = Bg + c * exp(-(X.^2 + Y.^2 + Z.^2 - (Z * ct1 + st1 * ...
    (X * ct2 + Y * st2)).^2) / (2 * sigmaPSF^2)) .* ( ...
    (erf(abs(tmp1) / c0) .* tmp1) ./ abs(tmp1) + ...
    (erf(abs(tmp2) / c0) .* tmp2) ./ abs(tmp2));
