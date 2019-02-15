function M = maxMatching(n,E)
% function M = maxMatching(n,E)
% maxMatching computes the maximum matching of a graph consisting of n
% nodes connected by the edges defined in E. A matching is a subset of the 
% edges defined in E. This matching maximizes the cardinality 
% (number of edges) of the subset so that each node is connected to another 
% node not more than once.
%
% Required Inputs:
% n              Scalar defining the number of nodes
% E              An m x 2 array representing m edges. Each line consists of 
%                two node indices. The order of the indices does not matter.
% 
% Optional Inputs:
%
% Outputs:
% M              A m x 1 logical array representing the matching. True
%                means the edge is part of the matching.
%
% Example: 
%                o - o
%                | X |
%                o - o
%
%                n = 4;
%                E = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
%                M = maxMatching(n,E)
%
%
% Sylvain Berlemont, Pascal Berard, October 2011