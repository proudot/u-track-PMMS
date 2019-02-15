function F = dLPoint2D(X, Y, xC, yC, A, Bg, sigmaPSF)
% Diffraction-limited Point Model
% F = dLPoint2D(X, Y, xC, yC, A, Bg, sigmaPSF);
%
% parameters:
% (X, Y)       2 vectors representing the 2-dimensional domain (e.g X =
%              -10:.1:10, Y = -5:.1:5
%
% (xC,yC)      center of the segment
%
% A            amplitude of the segment
%
% Bg           ampligure of the background
%
% sigmaPSF     half width of the gaussian PSF model.
%
% output:
% F is a NxM matrix where N = numel(X) and M = numel(Y).
%
% Sylvain Berlemont, 2009

[X Y] = meshgrid(X,Y);
X = X - xC;
Y = Y - yC;
F = Bg + A * exp(-0.5 * (X.^2 + Y.^2) / sigmaPSF^2);
end