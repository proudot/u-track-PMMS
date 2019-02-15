function C = pcombs2(v1,v2)
% This function generate all pairwise combinations between elements in v1
% with elements in v2. In both set, all elements must be different.
%
% e.g.:
% v1 = [1 2 3]
% v2 = [4 3]
% C = 
%   1  4
%   1  3
%   2  4
%   2  3
%   3  4
%   3  3

assert(length(v1) == length(unique(v1)) && length(v2) == length(unique(v2)));

[I J] = meshgrid(v1, v2');
C = [I(:), J(:)];

end