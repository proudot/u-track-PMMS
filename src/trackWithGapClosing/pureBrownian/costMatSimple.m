function [costMat,noLinkCost,nonlinkMarker,errFlag] = costMatSimple(...
    movieInfo,costMatParams)
%COSTMATSIMPLE provides a simple cost matrix for linking features between 2 time points
%
%SYNOPSIS [costMat,noLinkCost,nonlinkMarker,errFlag] = costMatSimple(...
%    movieInfo,costMatParams)
%
%INPUT  movieInfo    : A 2x1 array (corresponding to the 2 time points of 
%                      interest) containing the fields:
%             .xCoord      : Image coordinate system x-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%             .yCoord      : Image coordinate system y-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D. Default: zeros.
%             .zCoord      : Image coordinate system z-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D or 2D. Default: zeros.
%             .amp         : Amplitudes of PSFs fitting detected features. 
%                            1st column for values and 2nd column 
%                            for standard deviations.
%       costMatParams: Structure with the following fields:
%             .searchRadius: Maximum distance between two features in two
%                            consecutive time points that allows linking 
%                            them (in pixels).
%             .maxAmpRatio : Maximum ratio between the amplitudes of two
%                            features in two censecutive time points that 
%                            allows linking them.
%             .noLnkPrctl  : Percentile used to calculate the cost of
%                            linking a feature to nothing. Use -1 if you do
%                            not want to calculate this cost.
%
%OUTPUT costMat      : Cost matrix.
%       noLinkCost   : Cost of linking a feature to nothing, as derived
%                      from the distribution of costs.
%       nonlinkMarker: Value indicating that a link is not allowed.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The cost for linking feature i in time point t to feature j 
%in time point t+1 is given by
%(distance(i,j)^2)x(max(amp(i),amp(j))/min(amp(i),amp(j))).
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costMat = [];
noLinkCost = [];
nonlinkMarker = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatSimple')
    disp('--costMatSimple: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether problem is 1D, 2D or 3D and augment coordinates if necessary
if ~isfield(movieInfo,'yCoord') %if y-coordinates are not supplied

    %assign zeros to y and z coordinates
    for i=1:2
        movieInfo(i).yCoord = zeros(size(movieInfo(i).xCoord));
        movieInfo(i).zCoord = zeros(size(movieInfo(i).xCoord));
    end

else %if y-coordinates are supplied

    if ~isfield(movieInfo,'zCoord') %if z-coordinates are not supplied

        %assign zeros to z coordinates
        for i=1:2
            movieInfo(i).zCoord = zeros(size(movieInfo(i).xCoord));
        end

    end %(if ~isfield(movieInfo,'zCoord'))

end %(if ~isfield(movieInfo,'yCoord') ... else ...)

searchRadius = costMatParams.searchRadius;
maxAmpRatio = costMatParams.maxAmpRatio;
noLnkPrctl = costMatParams.noLnkPrctl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of features in the 2 time points
n = movieInfo(1).num;
m = movieInfo(2).num;

%replicate x,y-coordinates at the 2 time points to get n-by-m matrices
x1 = repmat(movieInfo(1).xCoord(:,1),1,m);
y1 = repmat(movieInfo(1).yCoord(:,1),1,m);
z1 = repmat(movieInfo(1).zCoord(:,1),1,m);
x2 = repmat(movieInfo(2).xCoord(:,1)',n,1);
y2 = repmat(movieInfo(2).yCoord(:,1)',n,1);
z2 = repmat(movieInfo(2).zCoord(:,1)',n,1);

%calculate the square distances between features in time points t and t+1
costMat = (x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2;

%assign NaN to all pairs that are separated by a distance > searchRadius
costMat(costMat>searchRadius^2) = NaN;

%replicate the feature amplitudes at the 2 time points to get n-by-m
%matrices
a1 = repmat(movieInfo(1).amp(:,1),1,m);
a2 = repmat(movieInfo(2).amp(:,1)',n,1);

%divide the larger of the two amplitudes by the smaller value
ampRatio = a1./a2;
for j=1:m
    for i=1:n
        if ampRatio(i,j) < 1
            ampRatio(i,j) = 1/ampRatio(i,j);
        end
    end
end

%assign NaN to all pairs whose amplitude ratio is larger than the
%maximum allowed
ampRatio(ampRatio>maxAmpRatio) = NaN;

%multiply the distance between pairs with the ratio between their
%amplitudes
costMat = costMat.*ampRatio;

%determine noLinkCost
if noLnkPrctl ~= -1
    noLinkCost = prctile(costMat(:),noLnkPrctl);
end

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;


%%%%% ~~ the end ~~ %%%%%

