function M = maxGreedyMatching(E,W)
% function M = maxGreedyMatching(E,W)
% maxGreedyMatching computes a greedy weighted matching of a graph 
% consisting of nodes connected by the edges defined in E and weighted 
% with W. A matching is a subset of the edges defined in E. Starting with 
% the edge with the highest weight, the edges, sorted by their corresponding 
% weights, are iteratively added to the matching as long as the edges do 
% not link to a node already linked by an edge already in the matching. 
% Each node is connected to another node not more than once. This function
% can be used as a greedy replacement function of the maximum matching functions.
%
% Required Inputs:
% E              An m x 2 array representing m edges. Each line consists of 
%                two node indices. The order of the indices does not matter.
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
%                M = maxWeightedMatching(E,W)
%
%
% Pascal Berard, December 2011

% Sort the edges
[~,sortIdx] = sort(W,'descend');
E = E(sortIdx,:);

M = false(size(W));
edgeIsActive = true(size(W));

while(any(edgeIsActive))
    
    % Find first active edge index
    tmp = find(edgeIsActive);
    firstActiveEdgeIdx = tmp(1);
    
    % Add the first active edge to the matching
    M(firstActiveEdgeIdx) = true;
    
    % Nodes of the first active edge
    firstActiveEdgeNodes = E(firstActiveEdgeIdx,:);
    
    % Deactivate all edges linking to a node of the first active edge
    deactivateActiveEdgeLeft = ismember(E(edgeIsActive,1),firstActiveEdgeNodes);
    deactivateActiveEdgeRight = ismember(E(edgeIsActive,2),firstActiveEdgeNodes);  
    deactivateActiveEdge = deactivateActiveEdgeLeft | deactivateActiveEdgeRight;
    edgeIsActive(edgeIsActive) = ~deactivateActiveEdge;
    
end

% "Un"-sort edges
M(sortIdx) = M;

end

