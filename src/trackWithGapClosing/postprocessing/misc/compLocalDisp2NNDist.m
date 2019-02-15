function [dispNNDistRatio] = compLocalDisp2NNDist(tracks)
%COMPLOCALDISP2NNDIST compares for every particle its displacement to its nearest neighbor distance
%
%SYNPOSIS [dispNNDistRatio] = compLocalDisp2NNDist(tracks)
%
%INPUT  tracks          : Output of trackCloseGapsKalman or trackWithGapClosing.
%
%OUTPUT dispNNDistRatio : Array with first column indicating particle
%                         displacement, second column indicating nearest
%                         neighbor distance, and third column indicating
%                         the ratio between the two. This information is
%                         collected on a frame by frame basis, and only
%                         particle pairs connected in consecutive frames
%                         contribute an entry.
%
%Khuloud Jaqaman, December 2009


%% Input

if nargin < 1
    disp('Please enter tracks')
    return
end

%convert tracks from structure to matrix format
if isstruct(tracks)
    tracksIn = tracks;
    tracks = convStruct2MatIgnoreMS(tracksIn);
end

%% calculation

%get the coordinates out of the track matrix
xCoord = tracks(:,1:8:end);
yCoord = tracks(:,2:8:end);
zCoord = tracks(:,3:8:end);

%calculate the displacement from one frame to the next
xDisp = diff(xCoord,1,2);
yDisp = diff(yCoord,1,2);
zDisp = diff(zCoord,1,2);
xyzDisp = sqrt( xDisp.^2 + yDisp.^2 + zDisp.^2 );

%get the nearest neighbor distances in each frame
nnDist = NaN(size(xyzDisp));
for iFrame = 1 : size(xyzDisp,2)
    xyzCoordCurrent = [xCoord(:,iFrame) yCoord(:,iFrame) zCoord(:,iFrame)];
    indxNonNaN = find(~isnan(xyzCoordCurrent(:,1)));
    if length(indxNonNaN) > 1
        xyzCoordCurrent = xyzCoordCurrent(indxNonNaN,:);
        distMat = createDistanceMatrix(xyzCoordCurrent,xyzCoordCurrent);
        distMat = sort(distMat,2);
        nnDist(indxNonNaN,iFrame) = distMat(:,2);
    end
end

%convert matrices into vectors and retain only non-NaN entries
xyzDisp = xyzDisp(:);
nnDist = nnDist(:);
indxNonNaN = find( ~isnan(xyzDisp) & ~isnan(nnDist) );
xyzDisp = xyzDisp(indxNonNaN);
nnDist = nnDist(indxNonNaN);
ratioDisp2Dist = xyzDisp ./ nnDist;

%% Output

dispNNDistRatio = [xyzDisp nnDist ratioDisp2Dist];


%% ~~~ the end ~~~
