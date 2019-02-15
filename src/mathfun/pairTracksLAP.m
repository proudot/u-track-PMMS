% links12 = pairTracksLAP(tracks1, tracks2)
% Links tracks from two different sets based on average distance between overlapping track segments.
%
% Inputs: tracks1, tracks2: Two structures with fields:
%
%         .x, .y, .start, .end
%

% Francois Aguet, September 16, 2010

function [links12 links21 dM] = pairTracksLAP(varargin)

if nargin==2
    tracks1 = varargin{1};
    tracks2 = varargin{2};
    
    n1 = length(tracks1);
    n2 = length(tracks2);
    
    % relative start/end positions   
    iStart1 = [tracks1.start];
    iEnd1 = [tracks1.end];
    iStart2 = [tracks2.start];
    iEnd2 = [tracks2.end];

    xVect1 = [tracks1.x];
    yVect1 = [tracks1.y];
    xVect2 = [tracks2.x];
    yVect2 = [tracks2.y];
    
elseif nargin==6   
    [xVect1, yVect1, iStart1, iEnd1, xVect2, yVect2, iStart2, iEnd2] = deal(varargin{:});
    n1 = numel(xVect1);
    n2 = numel(xVect2);
else
    error('Incompatible input arguments.');
end

lifetime1 = iEnd1 - iStart1 + 1;
lifetime2 = iEnd2 - iStart2 + 1;

pairIdx = pcombs2(1:n1, 1:n2);
% pidx1 = pairIdx(:,1);
% pidx2 = pairIdx(:,2);
% t1 = iStart1(pairIdx(:,1));
% t2 = iStart2(pairIdx(:,2));
% whos
% compute overlap between track pairs
overlapStart = max(iStart1(pairIdx(:,1)), iStart2(pairIdx(:,2)));
overlapEnd = min(iEnd1(pairIdx(:,1)), iEnd2(pairIdx(:,2)));
overlap = overlapEnd - overlapStart + 1;


% start, end indexes for these vectors
end1 = cumsum(lifetime1);
start1 = end1-lifetime1+1;
end2 = cumsum(lifetime2);
start2 = end2-lifetime2+1;


% translate overlap start/end values to 1-D indexes
start1 = start1(pairIdx(:,1)) + overlapStart - iStart1(pairIdx(:,1));
start2 = start2(pairIdx(:,2)) + overlapStart - iStart2(pairIdx(:,2));


% sort overlap values
[overlap idx] = sort(overlap);
pairIdx = pairIdx(idx,:);

start1 = start1(idx);
start2 = start2(idx);

minOverlap = 15;

startIdx = find([1 diff(overlap)] & overlap >= minOverlap);
endIdx = find([-diff(-overlap) 1] & overlap >= minOverlap);

% initialize cost matrix
dM = 100*ones(n1,n2);

overlapValues = unique(overlap(overlap>=minOverlap));


for k = 1:length(startIdx)
    
    % indexes corresponding to overlap value
    range = startIdx(k):endIdx(k);
    
    M = repmat((0:overlapValues(k)-1), [length(range) 1]);
    idx1 = repmat(start1(range)', [1 overlapValues(k)]) + M;
    idx2 = repmat(start2(range)', [1 overlapValues(k)]) + M;
    
    x1 = xVect1(idx1);
    y1 = yVect1(idx1);
    x2 = xVect2(idx2);
    y2 = yVect2(idx2);

    % average distance   
    meanDist = mean(reshape(sqrt((x1-x2).^2 + (y1-y2).^2), [length(range) overlapValues(k)]), 2);

    % fill cost matrix
    dM(sub2ind([n1 n2], pairIdx(range,1), pairIdx(range,2))) = meanDist;
end

nonLinkMarker = -1;
m12 = nonLinkMarker*ones(n1,n1);
maxSR = 5;
diagIdx = 1 + (0:n1-1)*(n1+1);
m12(diagIdx) = maxSR+1;

m21 = nonLinkMarker*ones(n2,n2);
diagIdx = 1 + (0:n2-1)*(n2+1);
m21(diagIdx) = maxSR+1;
%m22 = ones(n2,n1);

dM(dM>maxSR) = -1;

m22 = (max(dM(:))+1)*ones(n2,n1);


M = [dM m12; m21 m22];
[links12 links21] = lap(M, -1, [], 0, maxSR+1);
% [links12 links21] = lap(dM, -1, [], 1);

% usage with thresholded distance matrix
% dM(dM>5) = -1;
% whos dM
[links12x links21x] = lap(dM, -1, [], 1);

figure; plot(links12, 'k.'); hold on; plot(links12x, 'ro')