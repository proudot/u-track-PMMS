function [T R] = computeICP(X1, X2, numIter, tol) %#ok<STOUT,INUSD>
% [T R] = computeICP(X1, X2, numIter, tol)
%
% This functions compute the Iterative Closest Point algorithm (ICP).
% ICP finds the optimal translation and rotation between two free forms X1
% and X2, defined as set of 3-dimensional points. This function iteratively
% minimizes the following distance:
%
% dk = 1/n2 sum over {p_i in X2} ||x1_i - R * p_i - T||^2,
%
% where
% dk is the distance at iteration k,
% n2 = size(X2, 1),
% x1_i is the closest point to p_i belonging to X1.
%
% 'numIter' is maximum number of iterations for the Iterative Closest
% Point (ICP).
%
% 'tol' specify the precision under which dk - dk+1 is not significant
% anymore and the algorithm must stop.
%
% Ref: "A Method for Registration of 3-D Shapes", P. J. Besl and N. D.
% McKay, IEEE Transactions on Pattern Analysis and Machine Intelligence.
% Vol. 14, No. 2, 1992.
%
% Sylvain Berlemont, Dec 2009