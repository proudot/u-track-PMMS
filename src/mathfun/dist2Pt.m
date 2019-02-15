function [dist]=dist2Pt(coordsYX,ptYX)
%
% DESCRIPTION: Finds distance from multiple points to one point
%
% SYNOPSIS: [dist]=dist2Pt(coordsYX,ptYX)
%
% INPUT: 
%        coordsYX : (n x 2) or (n x 2 x s) matrix containing yx-coordinates
%                   
%        ptYX     : yx-coordinates of the point to which the distances are
%                   measured
%
% OUTPUT: 
%        dist     : matrix containing distance from each point to ptYX
%
% MATLAB VERSION (originally written on): 7.0.1.24704 (R14) Service Pack 1 Windows_NT
%
% USERNAME: kathomps
% DATE: 21-Apr-2006

% if coordsYX is n x 2, dist is n-vector of distances
% if coordsYX is n x 2 x s, dist is n x s matrix of distances

dist=squeeze(sqrt((coordsYX(:,1,:)-ptYX(1)).^2+(coordsYX(:,2,:)-ptYX(2)).^2));