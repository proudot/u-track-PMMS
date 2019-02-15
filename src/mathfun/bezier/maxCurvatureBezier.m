function [maxCurvature,tMax] = maxCurvatureBezier(cP,t0,t1)
% function [maxCurvature,tMax] = maxCurvatureBezier(cP,t0,t1)
% maxCurvatureBezier returns the value and the position of the maximum 
% curvature in the interval (t0,t1)
%
% Required Inputs:
% cP                A N x 3 array representing a set of 3-dimensional 
%                   control points where N is 2,3 or 4      
% 
% Optional Inputs:
% t0 = 0 (default)  Start value of the parametrization interval            
% t1 = 1 (default)  End value of the parametrization interval
%
% Outputs:
% maxCurvature      The maximum curvature
% tMax              The position of the maximum curvature
%
% Pascal Berard, June 2011