function fImg = locmin2d(img, mask, keepFlat)
%LOCMIN2D searches for local minima in an image
%
%    SYNOPSIS fImg = locmin2d(img, mask, keepFlat)
%
%    INPUT    img    image matrix
%             mask   EITHER a scalar that defines the window dimensions
%                    OR a vector [m n] that defines the window dimensions
%                    OR a binary (0/1) structural element (matrix).
%                    Structural elements such as discs can be defined
%                    using the built-in matlab function "strel".
%                    The input matrix must have an odd number of columns and rows. 
%             keepFlat Optional input variable to choose whether to remove
%                      "flat" maxima or to keep them. Default is 0, to remove them.
%
%    OUTPUT   fImg   image with local minima (original values) and zeros elsewhere.

if nargin < 3 || isempty(keepFlat)
    keepFlat = 0;
end

fImg = -locmax2d(-img, mask, keepFlat);