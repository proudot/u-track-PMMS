function maxGreedyMatchingTest

clear all; clc;

E = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
W = [2;6;9;8;50;3];

M = maxGreedyMatching(E,W)
E(M,:)
W(M)

end

