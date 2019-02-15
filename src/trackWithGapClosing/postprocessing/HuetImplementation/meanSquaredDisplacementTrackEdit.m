function mSqDisp = meanSquaredDisplacementTrackEdit(positions, sigmaZero, doOneD, analyzingParameters)
%MEANSQUAREDDISPLACEMENT calculates the mean squared displacement of coordinates that come with covariances
%
% SYNOPSIS mSqDisp = meanSquaredDisplacement(positions, sigmaZero, doOneD)
%
% INPUT    positions            1-by-(1-2) array of positions and covariances
%                               First column: position of tag.
%                               Second column (opt): position of reference
%               Fields:
%                   .coordinates    n-by-d array of coordinates
%                   .covariances    d-by-d-by-n array of covariances
%          
%          sigmaZero            (opt) n-by-1 array of sigma0 (noise) at the
%                               position of the coordinates. Default: 1
%
%          doOneD               (opt) whether or not to consider a
%                               higher-dimensional movement as 1d-motion
%                               [{0}/1]
%          
%          min
%               Please supply missing coordinates (and corresponding
%               covariances) as NaNs. Perfect observations have covariance
%               matrices of 0. (it will not turn out well if you mix perfect
%               and non-perfect observations!)
%
% OUTPUT   mSqDisp              x-by-3 array with
%                                   - r^2=<[r(t+dt)-r(t)]^2>
%                                   - sigma(r^2) - for SEM divide by 
%                                         sqrt(nDataPoints!)
%                                   - nDataPoints
%
%                               dt is max 1/4 of the
%                               available timepoints. 
%                               See Saxton, Biophys.J.,1997
%
%c: 05/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEFAULTS
defaultDoOneD = 0;
defaultSigmaZero = 1;
defaultFracTimepoints = 0.66;
defaultMinNumData = analyzingParameters.nMSD;

%=====================
% TEST INPUT
%=====================

% nargin 
if nargin < 1 || ~isstruct(positions) || ~isfield(positions,'coordinates') || ~isfield(positions,'covariances')
    error('Please specify a ''positions'' structure with fields ''coordinates'' and ''covariances''!')
end

% sizes
coordSize = size(positions(1).coordinates);
covSize   = size(positions(1).covariances);

% ntp
if coordSize(1) ~= covSize(3)
    error('Unequal number of coordinates and/or covariances')
end
% dims
if any(coordSize(2) ~= covSize(1:2))
    error('Dimension mismatch between coordinates and covariances')
end

% fill sigmaZero
if nargin < 2 | isempty(sigmaZero)
    sigmaZero = repmat(defaultSigmaZero,[coordSize(1),1]);
end

% fill referenceCoords, Covs if necessary
if length(positions) == 1 || isempty(positions(end).coordinates) || isempty(positions(end).coordinates)
    referenceCoordinates = zeros(coordSize);
    referenceCovariances = zeros(covSize);
else
    referenceCoordinates = positions(end).coordinates;
    referenceCovariances = positions(end).covariances;
end

% this is not so elegant, but I am quickly changing the code, so it will
% have to suffice. 
coordinates = positions(1).coordinates;
covariances = positions(1).covariances;
clear positions


if nargin < 3 || isempty(doOneD)
    doOneD = defaultDoOneD;
end

%========================


%========================
% DO 1-D IF SELECTED
%========================

if doOneD
    % calculate distance, sigma relative to reference, then use this for
    % new 1D coordinate/covariance matrix
    positions(1:2) = struct('coordinates',[],'covariances',[]);
    positions(1).coordinates = coordinates;
    positions(1).covariances = covariances;
    positions(2).coordinates = referenceCoordinates;
    positions(2).covariances = referenceCovariances;
    
    [distance,sigma] = deltaCoordinates(positions,sigmaZero);
    
    coordinates = distance;
    % make 1-by-1-by n array. Let reshape calculate the size
    covariances = reshape(sigma,1,1,[]);
    coordSize = size(coordinates);
    covSize   = size(covariances);
    referenceCoordinates = zeros(coordSize);
    referenceCovariances = zeros(covSize);
    
end

%=====================



%======================================================
% CALCULATE DISTANCES, FIND MEAN SQUARED DISPLACEMENTS
%======================================================

% find max time lag
maxTimeLag = floor(defaultFracTimepoints * length(find(isfinite(coordinates(:,1)))));

% preassign output
mSqDisp = repmat(NaN,[maxTimeLag,3]);

% calculate Input for deltaCoordinates. If all covariance matrices are
% zero, use unit matrices instead - we do not want to screw up our further
% calculations!
positions(1:3) = struct('coordinates',[],'covariances',[]);
positions(1).coordinates = coordinates;
if all(covariances == 0)
    covariances = repmat(eye(covSize(1)),[1,1,covSize(3)]);
end
positions(1).covariances = covariances;

positions(3).coordinates = referenceCoordinates;
if all(referenceCovariances == 0)
    referenceCovariances = repmat(eye(covSize(1)),[1,1,covSize(3)]);
end
positions(3).covariances = referenceCovariances;

% loop through all timeLags and calculate displacement and sigma
for timeLag = 1:maxTimeLag
    
    % calculate displacements.
    [displacement,sigma] = deltaCoordinates(positions,sigmaZero,timeLag);
    
%     % DEBUG
%     idx=find(sigma==0);
%     if ~isempty(idx)
%         idx,timeLag
%     end
    
    % find ~NaNs
    goodData = find(isfinite(displacement) & isfinite(sigma) & (sigma~=0));
    
    % if more than 20 data points, we average
    if length(goodData) >= defaultMinNumData
        % use sigma of distances, not squared distances, for weight
         [r2,dummy,sigmaR2] = weightedStats(displacement(goodData).^2,sigma(goodData),'s');
         mSqDisp(timeLag,:) = [r2,sigmaR2,length(goodData)];
    end
    
end % for timeLag = 1:maxTimeLag

