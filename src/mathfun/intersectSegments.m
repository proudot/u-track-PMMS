function [xi,yi] = intersectSegments(x1,y1,x2,y2)
%INTERSECTSEGMENTS finds the intersection of 2 line segments
%
% SYNOPSIS [xi,yi] = intersectSegments(x1,y1,x2,y2)
%
% INPUT x1 : vector with the x-coordinates of the 1st segment
%       y1 : vector with the y-coordinates of the 1st segment
%       x2 : vector with the x-coordinates of the 2nd segment
%       y2 : vector with the y-coordinates of the 2nd segment
%
% OUTPUT xi, yi : x/y coordinate pair of the intersection point
%
% NOTE  The intersection point can fall outside the segments.
%       The term segment is used as 2-point representation of a 
%       line from -infinity to +infinity


% check input data 
if(sum([length(x1),length(y1),length(x2),length(y2)]==2)~=4)
   error('invalid vectors entered');
end;

% build the equation system
A = [x1(2)-x1(1), x2(2)-x2(1);y1(2)-y1(1), y2(2)-y2(1)];
b = [x2(1)-x1(1);y2(1)-y1(1)];
s = A\b;

xi = x1(1) + s(1)*(x1(2)-x1(1));
yi = y1(1) + s(1)*(y1(2)-y1(1));
