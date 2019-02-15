function [dist,t] = distancePointBezier(cP,P,t0,t1)
% function [dist,t] = distancePointBezier(cP,P,t0,t1)
% distancePointBezier computes the distance between a point and a Bezier
% curve. If only one control point is specified the distance between the 
% point and the control point is returned. In comparison to distancePointBezierOld
% it is a bit slower. However, there is no limitation on the complexity of 
% the curves. This function makes use of the SISL NURBS library which can 
% be found at "http://www.sintef.no/sisl".
%
% Required Inputs:
% cP                A N x M array representing a set of N M-dimensional 
%                   control points.      
% 
% Optional Inputs:
% t0 = 0 (default)  Start value of the parametrization interval            
% t1 = 1 (default)  End value of the parametrization interval
%
% Outputs:
% dist              The distance to the Bezier curve
% t                 The corresponding curve parameter
%
% Pascal Berard, March 2012