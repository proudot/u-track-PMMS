function [fracClustTotalMSN,fracClustPerFrameMSN,fracClustTotalMovie,fracClustPerFrameMovie] = calcFracClust(tracksAggreg,maxSize)
%CALCFRACCLUST calculates the fraction of clusters of various sizes from series of tracks
%
%SYNOPSIS [fracClustTotalMSN,fracClustPerFrameMSN,fracClustTotalMovie,fracClustPerFrameMovie] = calcFracClust(tracksAggreg,maxSize)
%
%INPUT  tracksAggreg : Cell array of output of aggregStateFromCompTracks.
%                      1 entry = 1 movie or simulation.
%       maxSize      : Maximum cluster size to look into. Any clusters
%                      larger than maxSize will be lumped into the bin for
%                      maxSize.
%                      Optional. Default: maximum cluster size in data.
%
%OUTPUT fracClustTotalMSN     : maxSize-by-3 array with rows referring to
%                               clusters of size 1, 2, etc. and columns
%                               storing the mean over all frames over all
%                               movies, standard deviation over all frames
%                               over all movies, and number of movies
%                               contributing to the mean and std
%                               calculations.
%       fracClustPerFrameMSN  : maxSize-by-number of frames-by-3 array with
%                               rows referring to cluster size, columns
%                               referring to frames, and 3rd dimension
%                               storing, per frame, the mean over all
%                               movies, std over all movies, and number of
%                               movies.
%       fracClustTotalMovie   : Cell array with number of entries = number
%                               of movies. Each entry stores the cluster
%                               fractions over all frames in an individual
%                               movie.
%       fracClustPerFrameMovie: Cell array with number of entries = number
%                               of movies. Each entry is a
%                               maxSize-by-number of frames array stoaring
%                               the cluster fractions per frame in an
%                               individual movie.
%
%Khuloud Jaqaman, February 2013

%% Input

%get number of movies
numMovie = length(tracksAggreg);

%get matrices of aggregation state
aggregMat = cell(numMovie,1);
for iMovie = 1 : numMovie
    [~,~,~,~,aggregMat{iMovie}] = convStruct2MatIgnoreMS(tracksAggreg{iMovie}.defaultFormatTracks);
end

%get number of frames
aggregMatAll = vertcat(aggregMat{:});
numFrames = size(aggregMatAll,2);

%get maximum cluster size if not input
if nargin < 2 || isempty(maxSize)
    maxSize = max(aggregMatAll(:));
end

%% Cluster fractions

%fractions per frame in each movie
fracClustPerFrameMovie = cell(numMovie,1);
for iMovie = 1 : numMovie
    fracClust = hist(aggregMat{iMovie},1:maxSize);
    fracClust = fracClust ./ repmat(sum(fracClust,1),maxSize,1);
    fracClustPerFrameMovie{iMovie} = fracClust;
end

%overall fractions in each movie
fracClustTotalMovie = cell(numMovie,1);
for iMovie = 1 : numMovie
    fracClust = hist(aggregMat{iMovie}(:),1:maxSize);
    fracClust = fracClust / sum(fracClust);
    fracClustTotalMovie{iMovie} = fracClust;
end

%mean and std of fractions per frame
fracClustPerFrameAll = vertcat(fracClustPerFrameMovie{:});
fracClustPerFrameMSN = zeros(maxSize,numFrames,3);
for iSize = 1 : maxSize
    fracClustPerFrameMSN(iSize,:,1) = mean(fracClustPerFrameAll(iSize:maxSize:end,:));
    fracClustPerFrameMSN(iSize,:,2) = std(fracClustPerFrameAll(iSize:maxSize:end,:));
    fracClustPerFrameMSN(iSize,:,3) = numMovie;
end
    
%mean and std of overall fractions
fracClustTotalAll = vertcat(fracClustTotalMovie{:});
fracClustTotalMSN = [mean(fracClustTotalAll,1); std(fracClustTotalAll,[],1)]';
fracClustTotalMSN(:,3) = numMovie;

%% ~~~ the end ~~~
