function [gKernel] = anisotropicGaussD(sigmaX, sigmaY, theta, scaleFactor)
% Function "anisotropicGuassD" generates a 2D anisotropic Gauss filtering kernel, which can be used
% to convolve with the image to be filtered. The implmentation is straightforward and easy to read
% and change.
% 
% For better efficiency, use functions "anisotropicGaussC" or "anisotropicGaussR", 
% to be posted separately.
%  
% "anisotropicGaussC" implements the filter in a separable convolution kernel form.
% "anisotropicGaussR" implements the filter in a recursive filter form.
% SYNOPSIS [gKernel] = anisotropicGaussD(sigmaX, sigmaY, theta)
% 
% INPUT                  sigmaX  : variance along the X axis
%                        sigmaY  : variance along the Y axis
%                         theta  : in radian, rotation relative to the image coordinate system
%                   scaleFactor  : set the range of the kernel, e.g. set to 2 to assume 2 * sigma.
%                                  default value is 2.                         
% 
% OUTPUT                gKernel  : generated kernel of the anisotropic 2D filter
   
% DEPENDENCE        Requires only native matlab functions.
%   
% REMARKS           The filter kernel is not separated and therefore 2D. 
%   
% AUTHOR            Ge Yang, Ph.D.
%                   Laboratory for Computational Cell Biology
%                   The Scripps Research Institute
%                     
% DATE OF CREATION  March 08, 2004
%
% VERSION           0.2
%                   Minor changes to version 0.2 have been made.
%
% BUG REPORT        email: geyang@scripps.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Calculate the size of the 2D kernel

if nargin < 4 
    scaleFactor = 2;  % Assume 2 sigma by default. 
end
    
diagLength = scaleFactor * sqrt(sigmaX^2 + sigmaY^2);

% Calculate the size of the kernel
oriAngle = atan2(sigmaY, sigmaX);
newAngle1 = oriAngle + theta;
newAngle2 = -oriAngle + theta;
kernelWidthHalf = max(round(diagLength * abs(cos(newAngle1))), round(diagLength * abs(cos(newAngle2))));
kernelHeightHalf = max(round(diagLength * abs(sin(newAngle1))), round(diagLength * abs(sin(newAngle2))));

% Step 2: Compute filter coefficients
gKernel = zeros(2 * kernelHeightHalf + 1, 2 * kernelWidthHalf + 1);
for u = -kernelWidthHalf : 1 : kernelWidthHalf
    for v = -kernelHeightHalf : 1 : kernelHeightHalf
        % Compute the corresponding coordinate in the original (not
        % rotated) kernel
        x = u * cos(theta) + v * sin(theta);                    % u --> x
        y = u * sin(-theta) + v * cos(theta);                   % v --> y
        gKernel(v + kernelHeightHalf + 1, u + kernelWidthHalf + 1) = exp( - 0.5 * x^2 / (sigmaX^2) - 0.5 * y^2 / (sigmaY^2));
        % Since the coefficients will be normalized, it is not necessary to compute the coefficient of
        % 0.5/pi/sigma/sigma.
        
    end
end

% Step 3: Normalize the kernel
temp = sum(sum(gKernel));
gKernel = gKernel / temp;
