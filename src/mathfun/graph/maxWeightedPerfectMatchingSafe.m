function M = maxWeightedPerfectMatchingSafe(E,W)
% function M = maxWeightedPerfectMatchingSafe(E,W)
% maxWeightedPerfectMatching computes the maximum weighted perfect matching 
% of a graph consisting of n nodes connected by the edges defined in E and 
% weighted with W. A matching is a subset of the edges defined in E. This 
% matching maximizes the sum of the weights of the edges in the subset so 
% that each node is connected to another node EXACTLY ONCE. An even number 
% of nodes is required.
%
% The number of nodes as represented by the maximum value in E must be
% even.
%
% This safe version derives the number of nodes from the maximum value in
% E to prevent segmentation faults.
%
% Required Inputs:
% E              An m x 2 array representing m edges. Each line consists of 
%                two node indices. The order of the indices does not matter.
%                Each node needs to be connected at least with one other
%                node.
% W              An m x 1 array representing the weights of the m edges.
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
%                E = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
%                W = [2;6;7;8;5;3];
%                M = maxWeightedPerfectMatchingSafe(E,W)
%
%
% Mark Kittisopikul, Liya Ding, October 2014
M = maxWeightedPerfectMatching( max(E(:)),E,W);