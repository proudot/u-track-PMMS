function [dm,vectorMatrix] = distMat(M,metric)
%DISTMAT computes the distance matrix of a point set
%
% SYNOPSIS [dm,vectorMatrix] = distMat(M,metric);
%
% INPUT:  M: an nPoints x nDim matrix where the # of columns is the # of dimensions
%            and  the # of rows is the # of points
%         metric (optional): a nDim x nDim metric matrix (default: identity matrix)
%
% OUTPUT: dm: nPoints x nPoints, dm(2,1)=dm(1,2) is the distance between p1 and p2
%         vectorMatrix: nPoints x nPoints x nDim; vM(2,1,:) = -vM(1,2,:) is the
%                       vector from p2 to p1
%
% c: 20/04/01   dT

%get dims
[m, n] = size(M);

if nargin < 2
    metric=eye(n);
elseif length(metric)~= n
    error('incorrect dimension of metric');
end


if m < 2
   dm=0;
   return;
end

%calc indices so that only the necessary calculations are performed
p = (m-1):-1:2;
I = zeros(m*(m-1)/2,1);
I(cumsum([1 p])) = 1;
I = cumsum(I);
J = ones(m*(m-1)/2,1);
J(cumsum(p)+1) = 2-p;
J(1)=2;
J = cumsum(J);

%calc distances
Y = (M(I,:)-M(J,:))';
dm=squareform(sqrt(diag(Y'*metric*Y))');

if nargout == 2 %calculate vectors
    %prepare subscript matrix
    nVectors = length(J);
    [dimIdx,jiRowIdx,jiColIdx] = ndgrid(1:n,1:nVectors,1:2);
    jiRowIdx = jiRowIdx(:);
    jiColIdx = jiColIdx(:);
    dimIdx   = dimIdx(:);
    
    jiMatrix = [J,I];
    if size(jiMatrix,1) == 1
        jiMatrix = jiMatrix';
    end
    jiIdx_iM1 = sub2ind([nVectors,2],jiRowIdx,jiColIdx);
    jiIdx_iM2 = sub2ind([nVectors,2],jiRowIdx,jiColIdx(end:-1:1));
    
    %write matrix containing the indices into which the corresponding vector
    %entry is written (e.g. dimIdx = 1:nDim,1:nDim,1:nDim etc)
    indexMatrix = sub2ind([m,m,n],jiMatrix(jiIdx_iM1),jiMatrix(jiIdx_iM2),dimIdx);
    
    %init output matrix; for each entry i,j in dm, the corresponding vector is
    %vectorMatrix(i,j,:)
    vectorMatrix = zeros(m,m,n);
    
    %prepare vectors for reading
    vectors = [Y,-Y];
    vectors = vectors(:);
    
    %assign the vector components to vectorMatrix
    vectorMatrix(indexMatrix) = vectors;
    
else
    %no 2nd output argument, therefore do not do any calculations
    vectorMatrix = [];
    
end