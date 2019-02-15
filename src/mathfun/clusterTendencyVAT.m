function [sortedIndices,dmSort] = clusterTendencyVAT(dissimilarityMatrix)
%CLUSTERTENDENCYVAT visually assesses clustering tendency using the VAT algorithm
%
% SYNOPSIS: sortedIndices = clusterTendencyVAT(dissimilarityMatrix)
%
% INPUT dissimilarityMatrix : symmetric 2-d array of dissimilarities. It
%                             should have the following properties:
%                             dM(i,j) >= 0
%                             dM(i,j) == dM(j,i)
%                             dM(i,i) == 0
%
% OUTPUT sortedIndices: list of sorted rows/cols of the distance matrix
%        dmSort : sorted dissimilarity matrix. If dmSort is requested, the
%                 code will not display figures. 
%
% REMARKS Implementation of Bezdek & Hathaway 2002 VAT:A tool for visual
%         assessment of (Cluster) tendency
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: jdorn
% DATE: 12-Nov-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=======================
%% TEST INPUT
%=======================

if nargin == 0 || isempty(dissimilarityMatrix)
    error('Please input a non-empty dissimilarityMatrix')
end

% test for properties
if any(dissimilarityMatrix(:) < 0)
    error('No negative dissimilarities allowed')
end
if ~all(all(dissimilarityMatrix == dissimilarityMatrix',2),1)
    error('DissimilarityMatrix needs to be symmetric')
end
if ~all(diag(dissimilarityMatrix) == 0)
    error('Dissimilarity of objects with themselves have to be 0')
end

%======================



%=============================
%% ORDER DISSIMIARITYMATRIX
%=============================

% dm-size
nObjects = size(dissimilarityMatrix,1);

% nomenclature compared to article
% K: allIndices
% I: sortedList
% J: unsortedList
% P: sortedIndices
allIndices = (1:nObjects)';
sortedList = zeros(nObjects,1);
unsortedList = allIndices;
sortedIndices = zeros(nObjects,1);


% Find the largest dissimilarity
[i,j] = find(dissimilarityMatrix == max(dissimilarityMatrix(:)));
% there may be more than one value (which already gives two entries). Take
% the first one. Also, since the matrix is symmetric, it doesn't actually
% matter whether we sort rows or columns. Therefore, sort just rows
sortedIndices(1) = i(1);
sortedList(1) = i(1);
unsortedList(unsortedList == i(1)) = [];

% loop through nObjects-1 to find minimum spanning tree, which will give us
% the sortedIndices
for iRow = 2:nObjects
    % read all distances between grouped and ungrouped objects
    groupDistances = dissimilarityMatrix(sortedList(1:iRow-1),unsortedList);
    % find minimum distance between already sorted and the rest
    [i,j] = find(groupDistances == min(groupDistances(:)));
    % add to sortedIndices the jth unsorted column
    newRow = unsortedList(j(1));
    sortedIndices(iRow) = newRow;
    sortedList(iRow) = newRow;
    unsortedList(j(1)) = [];
end

%=================================

%=================================
%% DISPLAY UNSORTED, SORTED DM
%=================================

if nargout > 1
    dmSort = dissimilarityMatrix(sortedIndices,sortedIndices);
else

    fh = uiViewPanel;
    set(fh,'Name','Unsorted dissimilarity matrix');
    imshow(dissimilarityMatrix,[0,max(dissimilarityMatrix(:))]);

    fh = uiViewPanel;
    set(fh,'Name','Sorted dissimilarity matrix');
    imshow(dissimilarityMatrix(sortedIndices,sortedIndices),...
        [0,max(dissimilarityMatrix(:))]);
end