function [newX,newY] = rotate2D(theta,oldX,oldY)
%ROTATE2D : Compute the new coordinates after a 2D rotation about the origin.
%
% SYNOPSIS :
%    [newX newY] = rotate2D(theta,oldX,oldY)
%
% INPUT :
%    theta : The angle of rotation.
%    oldX  : The x-coordinates of the set of points to be rotated.
%    oldY  : The y-coordinates of the set of points to be rotated.
%
% OUTPUT :
%    newX  : The x-coordinates of the input points after rotation.
%    newY  : The y-coordinates of the input points after rotation.
%
% Lin Ji, Dec. 8, 2003

%Check inputs.
if ~isnumeric(theta) | length(theta) ~= 1
   error('The first input should be a scalar numerical value.');
end

if ~isnumeric(oldX) | ~isnumeric(oldY) | ...
   ndims(oldX) > 2 | ndims(oldY) > 2 | ...
   size(oldX,1) ~= size(oldY,1) | ...
   size(oldX,2) ~= size(oldY,2)
   error(['The 2nd and 3rd inputs must be numerical vectors of ' ...
      'the same length or 2D arrays of the same size.']);
end

%Construct the rotation matrix
R = [cos(theta) -sin(theta);sin(theta) cos(theta)];

[m n] = size(oldX);
oldX = reshape(oldX,1,m*n);
oldY = reshape(oldY,1,m*n);

newXY = R*[oldX;oldY];

newX = reshape(newXY(1,:),m,n);
newY = reshape(newXY(2,:),m,n);
