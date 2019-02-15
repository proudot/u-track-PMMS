function [q,g,h,r] = boundaryCondition(p,e,u,time)
%BOUNDARYGEOMETRY generic geometry M file for use with PDE toolbox.
%
% [q,g,h,r] = boundaryCondition(p,e,u,time)
%
% This function is designed for use with solvePDEVectorBoundary.m and the
% PDE toolbox.
% 
%
% Hunter Elliott
% Re-written 8/2010
%

%h & r specify dirichlet boundary conditions.

%Determine number of edges requested
nE = size(e,2);

global BOUND_COND;
global OBJ_BOUND;

%DIRICHLET CONSTRAINTS

%Function value multiplier.
h = [1 0 0 1]';
h = repmat(h,1,2*nE);

%Get parameter values for each edge requested         
s(1:nE) = e(3,1:nE);
s(1+nE:2*nE) = e(4,1:nE);
%Sometimes the pdetoolbox requests very small negative numbers (1e-20) we
%convert these to zero since the spline is undefined at negative parameter
%values
s(s<0) = 0;

%Now convert these to points spaced by arc length
iPars = linspace(0,1,OBJ_BOUND.pieces+1);
currEdge = fnval(OBJ_BOUND,iPars);
adjustedPar = pdearcl(iPars,currEdge,s,0,1);

%And evaluate the boundary condition at these points
r(:,1:nE) = ppval(BOUND_COND,adjustedPar(1:nE));
r(:,nE+1:2*nE) = ppval(BOUND_COND,adjustedPar(1+nE:2*nE));

%NEUMANN CONSTRAINTS
%Setting both to zero means dirichlet only

%function value multiplier
q = zeros(4,nE);

%Normal component of gradient at edge
g = zeros(2,nE);
