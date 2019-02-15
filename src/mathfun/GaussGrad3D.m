function [gaussGradY, gaussGradX, gaussGradZ] = GaussGrad3D(sigma,fSze,cent,cnorm,absCenter)
% GAUSSGRAD3D creates 3D Gaussian gradients
%
%    SYNOPSIS [gaussGradY, gaussGradX, gaussGradZ] = GaussGrad3D(sigma,fSze,cent,cnorm,absCenter)
%
%    INPUT: sigma  of gauss mask [sigmaX sigmaY sigmaZ] or [sigma]. 
%                  If scalar, sigma will be the same for all dimensions
%           fSze   size of the gauss mask [sizeX sizeY sizeZ] or [size]
%                  If scalar, size will be the same for all dimensions.
%                  (odd size required for symmetric mask!)
%           cent   (optional)3D vector with center position [0 0 0] is
%                  center of fSze (default=[0 0 0])
%                  if absCenter = 1, center position is measured from
%                  [0,0,0] in matrix coordinates (which is outside of the
%                  array by half a pixel!)
%           cnorm  (optional) select normalization method:
%                  =0 (default) no normalization - max of gauss will be 1
%                  =1 norm so that integral will be exactly 1
%                  =2 norm so that integral of an infinite Gauss would be 1
%           absCenter (optional) select the type of center vector
%                  =0 (default). Zero is at center of array
%                  =1 Zero is at [0,0,0] in matrix coordinates, which is
%                  outside of the array by half a pixel!
%                      (center of top left pixel of 2-D matrix = [1,1])
%
%
%    OUTPUT: gaussGradY, gaussGradX, gaussGradZ
%            gradients of 3D Gaussians
%            
%    REMARKS: In analogy to the built-in function gradient (which is slower
%               than convolution with central difference kernel),
%               GaussGrad3D returns first the derivation along the image
%               x-axis.
%
% 2/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
ls = length(sigma);
switch ls
    case 3
        % all is good
    case 1
        sigma = [sigma,sigma,sigma];
    otherwise
        error('sigma has to be either a 1-by-3 vector or a scalar!')
end

lf = length(fSze);
switch lf
    case 3
        % all is good
    case 1
        fSze = [fSze,fSze,fSze];
    otherwise
        error('fSze has to be either a 1-by-3 vector or a scalar!')
end
if nargin < 3 || isempty(cent)
    cent=[0 0 0];
end;
if nargin<4 || isempty(cnorm)
    cnorm=0;
end;
if nargin < 5 || isempty(absCenter)
    absCenter = 0;
end

% transform absolute center into relative center
if absCenter
    cent = cent-(fSze+1)/2;
end


% to get the accurate pixel intensities, use the cumsum of the gaussian
% (taken from normcdf.m). As the origin is in the center of the center
% pixel, the intensity of pixel +2 is the integral from 1.5:2.5, which,
% fortunately, has the spacing 1.

[gaussGradX,gaussGradY,gaussGradZ]=deal(zeros(fSze));

x=([-fSze(1)/2:fSze(1)/2]-cent(1))./sigma(1);
y=([-fSze(2)/2:fSze(2)/2]-cent(2))./sigma(2);
z=([-fSze(3)/2:fSze(3)/2]-cent(3))./sigma(3);
ex = diff(0.5 * erfc(-x./sqrt(2)));
ey = diff(0.5 * erfc(-y./sqrt(2)));
ez = diff(0.5 * erfc(-z./sqrt(2)));
% for gradient: use the same x, then use standard Gauss for differences
dex = diff(exp(-x.^2/2))/(sqrt(2*pi)*sigma(1));
dey = diff(exp(-y.^2/2))/(sqrt(2*pi)*sigma(2));
dez = diff(exp(-z.^2/2))/(sqrt(2*pi)*sigma(3));

% construct the 3D matrices
dxy = dex'*ey;
xdy = ex'*dey;
xy=ex'*ey;
gaussGradX(:) = -dxy(:) * ez;
gaussGradY(:) = -xdy(:) * ez;
gaussGradZ(:) = -xy(:) * dez;


% norm Gauss
switch cnorm
    case 0 % maximum of Gauss has to be 1 (unfortunately, grad(0)=0)
        gaussGradX = gaussGradX * ((2*pi)^1.5*prod(sigma));
        gaussGradY = gaussGradY * ((2*pi)^1.5*prod(sigma));
        gaussGradZ = gaussGradZ * ((2*pi)^1.5*prod(sigma));

    case 1
        gaussGradX = gaussGradX / sum(abs(gaussGradX(:)));
        gaussGradY = gaussGradY / sum(abs(gaussGradY(:)));
        gaussGradZ = gaussGradZ / sum(abs(gaussGradZ(:)));
        
    case 2 % the whole erfc thing is already normed for infinite Gauss
        % so nothing to do here.
        

    otherwise
        

end