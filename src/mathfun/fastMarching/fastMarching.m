function [U E S R] = fastMarching(X, F, maxDegree) %#ok<STOUT,INUSD>
% [U E S R] = FASTMARCHING(X, F, maxDegree) resolves the Eikonal Equation
% |grad(U)|F = 1, where U(X(:)) = 0. During this process, a graph G =
% (X, E) is created where an edge e is created between Xi and Xj as soon as
% their associated front, i.e. propagating from Xi and Xj, are touching
% each other for the first time.
%
% input:
%
%    X          initial set of point where fronts are propagating from. It
%               can be either a 2- or 3-dimensional point list. Coordinates
%               must be integers. Note that the convention for coordinate
%               order is X = <x, y> or <x, y, z>.
%
%    F          the speed function |grad(U)|F = 1. F must has the same
%               dimension (2- or 3-dimension) than X. Default value is
%               F(:) = 1 for every point.
%
%    maxDegree  it stands for the maximum number of edges for each vertex
%               in G.
%
% output:
%
%    U          the solution of the Eikonal Equation.
%
%    E          the set of edges between pair of points in X. It is a m x 2
%               matrix where X(E(i, 1), :) and X(E(i, 2), :) correspond to
%               the 2 points linked by the ith edge.
%
%    S          set of saddle points. A saddle point is the location where
%               2 fronts starting from Xi and Xj are touching for the first
%               time. The minimal path between Xi and Xj is therefore
%               obtained by concatenating the minimal path from S to Xi and
%               S to Xj. S is a m x 2 or m x 3 matrix (m is the number of
%               edges). S follows the same coordinate order convention,
%               i.e, <x, y> or <x, y, z>.
%
%    R          the front map. R(i, j) = n means that the nth front
%               propagating from X(n, :) is the first front reaching the
%               pixel (i, j). Same definition in 3D.
%
% Reference:
% "Level Set Methods and Fast Marching Methods: Evolving Interfaces in
% Computational Geometry, Fluid Mechanics, Computer Vision and Materials
% Science" J.A. Sethian, Cambridge University Press, 1999.
%
% Sylvain Berlemont, 2009
