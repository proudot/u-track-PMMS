function [featDensity,errFlag] = calcFeatDensityPerMotionType(trackData,...
    fullImagePathName,diffAnalysisData,minTrackLen)


%% input

errFlag = 0;

if nargin < 3
    disp('calcFeatDensityPerMotionType: Missing input argument!');
    return
end

if nargin < 4 || isempty(minTrackLen)
    minTrackLen = 5;
end

%% calculation

%get number of input movies
numMovies = length(trackData);

%reserve memory for output
featDensity = NaN(numMovies,6);

%go over all movies in trackData ...
for iMovie = 1 : numMovies
    
    %% preamble
    
    %get this movie's tracks and their diffusion analysis results
    tracks = trackData(iMovie).tracks;
    diffAnalysisRes = diffAnalysisData(iMovie).diffAnalysisRes;
    
    %keep only tracks with length >= minTrackLen
    criteria.lifeTime.min = minTrackLen;
    indx = chooseTracks(tracks,criteria);
    clear criteria
    tracks = tracks(indx);
    diffAnalysisRes = diffAnalysisRes(indx);
    
    %put tracks in matrix format
    [tracksMat,tracksIndxMat] = convStruct2MatIgnoreMS(tracks);
    
    %get number of frames
    numFrames = size(tracksIndxMat,2);
    
    %load image
    cellImage = imread(fullImagePathName{iMovie});
    
    %get area of image (i.e. cell area) used for feature detection
    cellArea = length(find(cellImage~=0));
    
    %% motion types
    
    %get track segment classification from diffusion analysis
    trackSegmentClass = vertcat(diffAnalysisRes.classification);
    
    %get track segment length
    trackSegmentLength = getTrackSEL(tracksMat);
    trackSegmentLength = trackSegmentLength(:,3);
    
    %get indices of linear, Brownian, confined Brownian and undetermined tracks
    indxLin    = find( trackSegmentClass(:,1) == 1 | trackSegmentClass(:,2) == 3 );
    indxBrown  = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 2 );
    indxConf   = find( trackSegmentClass(:,1) ~= 1 & trackSegmentClass(:,2) == 1 );
    indxUndet1 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
        & trackSegmentLength >= 5);
    indxUndet2 = find( trackSegmentClass(:,1) ~= 1 & isnan(trackSegmentClass(:,2)) ...
        & trackSegmentLength < 5);
    
    %calculate average number of features in each category per frame
    numFeatLin    = length(find(tracksIndxMat(indxLin,:))) / numFrames; %#ok<FNDSB>
    numFeatBrown  = length(find(tracksIndxMat(indxBrown,:))) / numFrames; %#ok<FNDSB>
    numFeatConf   = length(find(tracksIndxMat(indxConf,:))) / numFrames; %#ok<FNDSB>
    numFeatUndet1 = length(find(tracksIndxMat(indxUndet1,:))) / numFrames; %#ok<FNDSB>
    numFeatUndet2 = length(find(tracksIndxMat(indxUndet2,:))) / numFrames; %#ok<FNDSB>
    
    %get density of features in each category
    densityLin    = numFeatLin    / cellArea;
    densityBrown  = numFeatBrown  / cellArea;
    densityConf   = numFeatConf   / cellArea;
    densityUndet1 = numFeatUndet1 / cellArea;
    densityUndet2 = numFeatUndet2 / cellArea;
    
    %put all densities in one row
    densityAll = [densityLin densityBrown densityConf densityUndet1 densityUndet2];
    
    %calculate overall density
    densityAll = [sum(densityAll) densityAll];
    
    %put this movies densities in their proper row
    featDensity(iMovie,:) = densityAll;
    
end

