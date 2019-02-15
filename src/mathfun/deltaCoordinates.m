function [distance, sigma, unitVector, sigmaUnitVector] = deltaCoordinates(points,sigmaZero,timeLag)
%DELTACOORDINATES calculates distance, uncertainty and unit orientation vector from points with corresponding uncertainty
%
% SYNOPSIS [distance, sigma, unitVector] = deltaCoordinates(points,sigmaZero,timeLag)
%
% INPUT    points    structure array of length 1-4 with coordinates and covariances
%                       of the individual points. n is the amount of e.g.
%                       timepoints for which you want to measure the
%                       distances. The second dimension is ordered
%                       [point1, point2, correction1, correction2]. Points
%                       are the point coordinates, corrections are the
%                       coordinates by which the points have to be shifted.
%                       If you want to calculate a displacement, fill
%                       column 1 (and 3 for corrected). If you want
%                       to calculate a distance measurement between two
%                       points for all n timepoints, fill col. 1 and 2. If
%                       both points have to be shifted by different
%                       centers, fill col. 1 to 4. You will always need to
%                       fill col. 1. Fill in zeros if e.g. a reference
%                       point is perfectly known
%             fields of "points" (n: number of timepoints, d: dimension)
%               .coordinates : nxd coordinates
%               .covariances : dxdxn array of covariances (or nxd array of
%                              diagonals of covariance matrices)
%             Please specify missing data points with NaN-entries of the
%               correct size
%
%          sigmaZero  (opt) n-by-1 array of average sigma0's (noise) of the
%                       measurements. Default: 1
%
%          timeLag    (opt) if you want to calculate displacements for a
%                       lag > 1, specify it here. Default: 1. Input for
%                       timeLag is not considered if there is data in the
%                       second col. of points.
%
%
% OUTPUT   distance    (n-timeLag+1)-by-d array of distances/displacements
%
%          sigma       uncertainty of the distances
%
%          unitVector  unit vector pointing from the first to the second
%                      point
%           
%          sigmaUnitVector uncertainty of every single component of the
%                          unit vector (normed to vector of length 1!)
%
%c: jonas 05/04 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEFAULTS
defaultSigmaZero = 1;
defaultTimeLag   = 1;

%=====================
% TEST INPUT
%=====================

% test nargin, isstruct
if nargin == 0 | ~isstruct(points) | ~isfield(points,'coordinates') | ~isfield(points,'covariances')
    error('DELTACOORDINATES needs a structure with fields ''coordinates'' and ''covariances'' as input')
end




coordSize = size(points(1).coordinates);
covSize   = size(points(1).covariances);

% assign sigmaZero, timeLag
if nargin < 2 | isempty(sigmaZero)
    sigmaZero = repmat(defaultSigmaZero,[coordSize(1),1]);
end
if nargin < 3 | isempty(timeLag)
    timeLag = defaultTimeLag;
end
if any([coordSize,covSize]==0)
    error('please specify at least coordinates for the first point')
end

% test dimension
if coordSize(2) ~= covSize(1) | coordSize(2) ~= covSize(2)
    if all(coordSize == covSize)
        % only the diagonal elements of the covariance matrices are
        % supplied. Make covariance matrices
        for iPoint = 1:length(points)
            if ~isempty(points(iPoint).covariances)
                covSizeOld = size(points(iPoint).covariances);
                covNew = zeros(covSizeOld(2),covSizeOld(2),covSizeOld(1));
                for t = 1:covSizeOld(1)
                    covNew(:,:,t) = diag(points(iPoint).covariances(t,:));
                end
                points(iPoint).covariances = covNew;
            end
        end

        % measure covSize again
        covSize = size(points(1).covariances);

    else
        error('coordinate dimensions do not match dimensions of covariances!')
    end
end
% test ntp - account for single distance
if ~(coordSize(1) == 1 && length(covSize) == 2) && coordSize(1) ~= covSize(3)
    error('number of coordinates does not match number of covariances!')
end

%---- assign distanceVector and all the covariances. If they are zeros,
%     they will not result in anything in the subsequent calculations
[distanceVectors,cov1,cov2,cov3,cov4] = fillProgramVars(points,sigmaZero,timeLag,coordSize,covSize);


%===================================================



%==========================
% CALC DISTANCES
%==========================

% call normList
[distance,unitVector] = normList(distanceVectors);

%=========================


%=========================
% CALC SIGMAS
%=========================

% because of matrix multiplications, we have to loop.
% What we will be doing: sigmaD = sqrt(sigma0*Qdd), Qdd = HQH',
% H=[-e,e,e,-e], Q=blkdiag(cov1,cov2,cov3,cov4)

if nargout > 1 % speed up stuff

    % assign Q
    nDims = coordSize(2);
    Q = zeros(4*nDims);

    % assign sigma
    nSigmas = size(distance,1);
    sigma = zeros(nSigmas,1);


    % loop
    for iSigma = 1:nSigmas

        % Q
        Q(0*nDims+1:1*nDims, 0*nDims+1:1*nDims) = cov1(:,:,iSigma);
        Q(1*nDims+1:2*nDims, 1*nDims+1:2*nDims) = cov2(:,:,iSigma);
        Q(2*nDims+1:3*nDims, 2*nDims+1:3*nDims) = cov3(:,:,iSigma);
        Q(3*nDims+1:4*nDims, 3*nDims+1:4*nDims) = cov4(:,:,iSigma);

        % H
        H = [-unitVector(iSigma,:), unitVector(iSigma,:), unitVector(iSigma,:), -unitVector(iSigma,:)];

        % sigma
        sigma(iSigma) = sqrt(sigmaZero(iSigma) * (H*Q*H'));
        
        if nargout > 3
            % calculate the uncertainties in the individual directions
            sigmaUnitVector = zeros(nSigmas,nDims);
            for d=1:nDims
                factor = [-1,1,1,-1];
                H = zeros(1,4*nDims);
                H(d:nDims:end) = factor;
                % divide sigma by the length of the distanceVector to get
                % the sigma of the components of the unit vector
                sigmaUnitVector(iSigma,d) = ...
                    sqrt(sigmaZero(iSigma) * (H*Q*H')) / distance(iSigma);
            end
        end

    end % for iSigma = 1:nSigmas

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% subfunctions
function [distanceVectors,cov1,cov2,cov3,cov4] = fillProgramVars(points,sigmaZero,timeLag,coordSize,covSize)
% fills in all the variables needed for computation



% we now need to fill 1 distance vector array and 4 covariance arrays. Select
% which ones and how
which2fill = (length(points)>1 && ~isempty(points(2).coordinates))+ ...
    2*(length(points)>2 && ~isempty(points(3).coordinates)) + ...
    4*(length(points)>3 && ~isempty(points(4).coordinates));

switch which2fill
    % cases
    % 7: all are full. shift pos1 and pos2, retain covariances
    % 3/5: only pos1/pos2 is shifted. Currently not supported (but easy to implement)
    % 2: displacement. Shorten pos1 by one, fill pos2 with t+1, correct
    %    both
    % 1: distance, no correction. Do as in 7, don't correct
    % 0: displacement, no correction. Do as in case 2, don't correct
    % 4/6: error

    case 7
        % shift both positions, retain all covs

        if all(size(points(2).coordinates)==coordSize) & all(size(points(3).coordinates)==coordSize) & all(size(points(3).coordinates)==coordSize) &...
                all(size(points(2).covariances)==covSize) & all(size(points(3).covariances)==covSize) & all(size(points(4).covariances)==covSize)

            distanceVectors = points(2).coordinates-points(4).coordinates - (points(1).coordinates-points(3).coordinates);
            cov1 = points(1).covariances;
            cov2 = points(2).covariances;
            cov3 = points(3).covariances;
            cov4 = points(4).covariances;

        else
            error('Unequal numbers of coordinates or covariances')
        end

    case 1
        % distance, no correction

        if all(size(points(2).coordinates)==coordSize) &...
                all(size(points(2).covariances)==covSize)

            distanceVectors = points(2).coordinates - points(1).coordinates;
            cov1 = points(1).covariances;
            cov2 = points(2).covariances;
            cov3 = zeros(covSize);
            cov4 = zeros(covSize);

        else
            error('Unequal numbers of coordinates or covariances')
        end

    case 2
        % displacement, corrected.

        if all(size(points(3).coordinates)==coordSize)  &...
                all(size(points(3).covariances)==covSize)

            pos1tmp = points(1).coordinates-points(3).coordinates;
            pos1 = pos1tmp(1:end-timeLag,:);
            pos2 = pos1tmp(timeLag+1:end,:);
            distanceVectors = pos2-pos1;
            cov1 = points(1).covariances(:,:,1:end-timeLag);
            cov2 = points(1).covariances(:,:,timeLag+1:end);
            cov3 = points(3).covariances(:,:,1:end-timeLag);
            cov4 = points(3).covariances(:,:,timeLag+1:end);

            % shorten sigmaZero
            sigmaZero = mean(sigmaZero(1:end-timeLag)+sigmaZero(timeLag+1:end),2);

            clear pos1tmp pos1 pos2

        else
            error('Unequal numbers of coordinates or covariances')
        end

    case 0
        % displacement, uncorrected

        pos1tmp = points(1).coordinates;
        pos1 = pos1tmp(1:end-timeLag,:);
        pos2 = pos1tmp(timeLag+1:end,:);
        distanceVectors = pos2-pos1;
        cov1 = points(1).covariances(:,:,1:end-timeLag);
        cov2 = points(1).covariances(:,:,timeLag+1:end);
        cov3 = zeros(covSize);
        cov4 = zeros(covSize);

        % shorten sigmaZero
        sigmaZero = mean(sigmaZero(1:end-timeLag)+sigmaZero(timeLag+1:end),2);

        clear pos1tmp pos1 pos2

    otherwise

        error('Please fill either cols 1, 1&2, 1&3 or 1&2&3&4 in points-structure')
end