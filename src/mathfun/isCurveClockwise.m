function isClockwise = isCurveClockwise(curveIn)
%ISCURVECLOCKWISE determines whether the points in the input curve run clockwise or counter-clockwise
%
% isClockwise = isCurveClockwise(curveIn)
%
% Returns true if the points in the closed, 2-D input curve run clockwise
% or counter clockwise in a right-handed coordinate system. (Will be
% reversed from matrix/image coordinates!)
%
% Input:
% 
%   curveIn - 2xM or Mx2 curve to check the orientation of. 
% 
% Output:
%
%   isClockwise - True if the points run clockwise, and false otherwise. 
% 
% 
%Hunter Elliott
%2008

if ndims(curveIn) ~= 2
    error('input curve of wrong dimension!');    
end

if length(curveIn) < 3
    error('Need at least three points to make a closed curve!')    
end

%Transpose to row vector if necessary
if size(curveIn,1) > size(curveIn,2)
    curveIn = curveIn';
end

%use the polygon area formula to check curve orientation
v1 = curveIn(1,[end 1:end-1]) .* curveIn(2,[2:end 1]);
v2 = curveIn(1,[2:end 1]) .* curveIn(2,[end 1:end-1]);

%If the area is negative, the curve is clockwise.
isClockwise = (sum(v1 - v2) / 2) < 0 ;