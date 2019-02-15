function [pts,values] = gradientDescent(F,X,Y) %#ok<STOUT,INUSD>
% Usage: [pts,values] = gradientDescent(F,X,Y)
%
% This function returns the paths along the steepest gradient of a 2D
% function F, starting from point X,Y. The number of paths returns is
% equivalent to size(X,1). vectors X and Y must be the same size.
%
% Outputs:
% pts          a cell array of size(X,1) elements. Each element is a 2D
%              [Mi X 2] matrix representing the list of Mi points along the
%              path.
%
% values       a cell array of size(X,1) elements where each element is a
%              vector of Mi elements representing the subpixellic value of
%              function F. Subpixellic values are obtained by bilinear
%              interpolation.
%
% e.g:
%
% % Get a mask
% mask = imread(...);
% % Compute the distance transfrom from the mask outline
% D = double(bwdist(1 - double(mask)));
% % Extract the mask outline
% C = contourc(D, [0, 0]);
% % Get a set of points lying on C
% x = C(1,5:5:end);
% y = C(2,5:5:end);
% % Get the steepest gradient descent paths
% [pts,values] = gradientDescent(max(D(:))-D,x,y);
% % Display
% imshow(D, []);
% hold on;
% n = c(2, 1);
% line(c(1, 2:n+1), c(2, 2:n+1), 'Color', 'y');
% for i = 1:numel(pts)
% line(pts{i}(:,1), pts{i}(:,2), 'Color', 'r');
% end