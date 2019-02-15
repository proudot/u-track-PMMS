function [trackClass,diffCoef,confParam,asymParam] = trackMSDAnalysis(...
    tracks,probDim)
%TRACKMSDANALYSIS classifies trajectories as pure Brownian, confined Brownian and Brownian with linear drift based on their mean square displacement
%
%SYNPOSIS [trackClass,diffCoef,confParam,asymParam] = trackMSDAnalysis(...
%    tracks,probDim)
%
%INPUT  tracks    : Matrix indicating the positions and amplitudes of the
%                   tracked features. Number of rows = number of tracks,
%                   number of columns = 8*number of time points. 
%                   Each row consists of 
%                   [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...].
%                   NaN is used to indicate time points where the track
%                   does not exist.
%       probDim   : Problem dimensionality. Optional. Default: 2.
%
%OUTPUT trackClass: # tracks x 1 vector of track classification.
%                   Values mean the following ...
%                   0 = stalled
%                   1 = confined Brownian
%                   2 = pure Brownian
%                   3 = Brownian with drift (directed)
%                   NaN = not classified
%       diffCoef  : # tracks x 1 vector of diffusion coefficients. NaN
%                   indicates tracks that could not be analyzed.
%       confParam : # tracks x 1 vector of parameter quantifying deviation
%                   of MSD curve from the straight line of a Brownian
%                   process. Values < 0 indicate confined motion (used for
%                   classification). Values > 0 indicate directed motion
%                   (not used for classification).
%                   NaN indicates tracks that could not be analyzed.
%       asymParam : # tracks x 1 vector of parameter quantifying asymmetry
%                   in particle positions along track. Larger values
%                   indicate higher asymmetry hence directed motion (used
%                   for classification).
%                   NaN indicates tracks that could not be analyzed.
%
%REMARKS (1) Algorithm is based on Huet at al. 2006. Biophys J 91:
%3542-3559. (2) No time-windowing for now, i.e. only whole trajectories are
%classified.
%
%Khuloud Jaqaman, February 2008

%% input

if nargin < 1 
    disp('--trackMSDAnalysis: Please input tracks to analyze.');
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

%get number of tracks and frames
numTracks = size(tracks,1);

%find indices of tracks that are >= 20 frames long - do not attempt
%to get diffusion coefficients for shorter tracks
criteria.lifeTime.min = 20;
indx4diff = chooseTracks(tracks,criteria);
clear criteria

%find indices of tracks whose length >= 5 and < 20 
%for asymmetry calculation, 5 frames and above are OK
% criteria.lifeTime.min = 5;
criteria.lifeTime.min = 5;
criteria.lifeTime.max = 19;
indx4asym = chooseTracks(tracks,criteria);
clear criteria

%% memory for trajectory classification

%classification means ...
%0 = stalled
%1 = confined Brownian
%2 = pure Brownian
%3 = drift/directed
%NaN = unclassified

trackClass = NaN(numTracks,1);
diffCoef = NaN(numTracks,1);
confParam = NaN(numTracks,1);
asymParam = NaN(numTracks,1);

%% confined Brownian motion

%define number of time points to use in calculating diffusion coefficient
nPointsCoef = 5;

%assign alpha-value for threshold determination
alphaConf = 0.1;

%go over all trajectories where in principle we should be able to calculate
%a mean square displacement curve and hence a diffusion coefficient and
%deviation parameter
for iTrack = indx4diff'
   
    %determine whether the track's MSD deviates (negatively) from that of a
    %pure Brownian track
    [confParamT,confFlag,diffCoefT] = msdDeviation2D3D(tracks(iTrack,:),...
        nPointsCoef,alphaConf,probDim);

    %classify track as ...
    %1 = confined Brownian, if the deviation is negative and its magnitude
    %is larger than the threshold
    %otherwise, keep track classification as undetermined
    if confFlag == 1
        trackClass(iTrack) = 1;
    end
    
    %save diffusion coefficient and deviation parameter
    diffCoef(iTrack) = diffCoefT;
    confParam(iTrack) = confParamT;
    
end

%find indices of tracks not classified yet
indxNotClass = find(isnan(trackClass));

%find indices of tracks where diffusion coefficient could not be calculated
indxNoDiffCoef = intersect(indx4diff,find(isnan(diffCoef)));

%remove these tracks from the list indx4diff and add them to the list
%indx4asym
indx4diff = setdiff(indx4diff,indxNoDiffCoef);
indx4asym = [indx4asym; indxNoDiffCoef];

%% directed motion

%assign alpha-value for threshold determination
alphaAsym = 0.1;

%go over all trajectories in indx4diff that haven't been classified yet
for iTrack = intersect(indx4diff',indxNotClass)

    %get the particle positions along the track
    coordX = tracks(iTrack,1:8:end)';
    coordY = tracks(iTrack,2:8:end)';
    coordZ = tracks(iTrack,3:8:end)';
    coordXYZ = [coordX coordY coordZ];
    
    %determine whether the track is sufficiently asymmetric
    [asymParamT,asymFlag] = asymDeterm2D3D(coordXYZ(:,1:probDim),alphaAsym);

    %classify track as ...
    %3 = directed, if the asymmetry parameter is larger than the threshold
    %2 = pure Brownian, if the asymmetry parameter is smaller than the
    %threshold (remember this track's deviation parameter was also smaller
    %than the deviation threshold)
    %otherwise, keep track classification as undetermined
    trackClass(iTrack) = 2 + asymFlag;

    %save asymmetry parameter
    asymParam(iTrack) = asymParamT;
    
end

%also go over trajectories that were too short to calculate a diffusion
%coefficient for them
for iTrack = indx4asym'

    %get the particle positions along the track
    coordX = tracks(iTrack,1:8:end)';
    coordY = tracks(iTrack,2:8:end)';
    coordZ = tracks(iTrack,3:8:end)';
    coordXYZ = [coordX coordY coordZ];
    
    %determine whether the track is sufficiently asymmetric
    [asymParamT,asymFlag] = asymDeterm2D3D(coordXYZ(:,1:probDim),alphaAsym);

    %classify track as ...
    %3 = directed, if the asymmetry parameter is larger than the threshold
    %otherwise, keep track classification as undetermined
    if asymFlag == 1
        trackClass(iTrack) = 3;
    end

    %save asymmetry parameter
    asymParam(iTrack) = asymParamT;
    
end

%% output


%% ~~~ the end ~~~