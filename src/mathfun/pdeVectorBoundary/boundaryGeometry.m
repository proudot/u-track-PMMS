function [x,y] = boundaryGeometry(bs,s)
%BOUNDARYGEOMETRY generic geometry M file for use with PDE toolbox.
%
% [x,y] = boundaryGeometry(bs,s)
%
% This function is designed for use with solvePDEVectorBoundary.m and the
% PDE toolbox.
% 
%
% Hunter Elliott
% Re-written 8/2010
%

%We need the boundary coordinates to be a global variable, because the PDE
%toolbox won't pass any extra arguments and won't accept a function handle
%as an input
global OBJ_BOUND;

%One segment per edge segment to avoid problems with strange geometries.
nBoundSeg = OBJ_BOUND.pieces;

%If no arguments supplied, return the number of boundary segments. This is
%standard for geometry m file
if nargin == 0
    x = nBoundSeg;
    return
end

%Since we expect input curve to be clockwise, this will label inside 1.
dirVec = [0 1];

%Initialze d matrix which defines the boundary segments
d = nan(4,nBoundSeg);

for j = 1:nBoundSeg
    %The parameter runs from 0 to 1 along the cell border.
    %The interior is labelled 1 and the exterior zero
    d(:,j) = [ (j-1)/nBoundSeg j/nBoundSeg dirVec ]';        
end

%Standard output if 1 argument for geometry m-file
if nargin == 1
    x = d(:,bs);
    return
end

%Adjust so points are evenly spaced by arc length. Only necessary 
%if there are gaps in boundary, remove if speed needed
iPars = linspace(0,1,OBJ_BOUND.pieces+1);
currEdge = fnval(OBJ_BOUND,iPars);
adjustedPar = pdearcl(iPars,currEdge,s,0,1);

% Evaluate this spline at the requested points 
%yy = ppval(OBJ_BOUND,s);
yy = ppval(OBJ_BOUND,adjustedPar);

x = nan(size(s));
y = nan(size(s));

%Split up the x & y coord for return
x(:) = yy(1,:);
y(:) = yy(2,:);


