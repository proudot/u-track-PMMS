function [movieInfo,tracksSimMiss,percentMissingActual] = genMovieInfoFromTracks2(tracksSim,...
    percentMissing,percentFP,minSegLength)
%GENMOVIEINFOFROMTRACKS2 generates a list of detected features per frame from supplied tracks
%
%SYNOPSIS [movieInfo,tracksSimMiss] = genMovieInfoFromTracks2(tracksSim,percentMissing,percentFP)
%
%INPUT  tracksSim     : Output of simulateMimickCD36_MS.
%       percentMissing: Percentage of missing features in movie.
%       percentFP     : Percentage of false detecton positives, relative to
%                       original number of features.
%                       Optional. Default: 0.
%       minSegLength  : Minimum length of a track segment between two gaps.
%                       Optional. Default: 1.
%
%OUTPUT movieInfo: List of detected features per frame, in the format
%                  required for the input of trackWithGapClosing and
%                  trackCloseGapsKalman.
%
%Khuloud Jaqaman, October 2007

if nargin < 3 || isempty(percentFP)
    percentFP = 0;
end

if nargin < 4 || isempty(minSegLength)
    minSegLength = 1;
end

%get number of frames in movie
seqOfEvents = vertcat(tracksSim.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%get number of tracks
numTracks = length(tracksSim);

%define standard deviation of missing features
missStd = 0.1*percentMissing;

%define standard deviation of false positives
fpStd = 0.1*percentFP;

%pre-allocate memory for movieInfo
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numFrames,1);

%get number of segments making each track
numSegments = zeros(numTracks,1);
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksSim(iTrack).tracksCoordAmpCG,1);
end

%locate the row of the first segment of each compound track in the
%big matrix of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%% put all tracks together in a big matrix
%store coordinates just before/after a merge/split as their negative values
%in order to distinguish them from coordinates not related to merges/splits
%also store in negative values the coordinates at the start/end of a
%segment
trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numFrames);

for iTrack = 1 : numTracks
    
    %get track's sequence of events and coordinates
    seqOfEvents = tracksSim(iTrack).seqOfEvents;
    tracksCoordAmpCG = tracksSim(iTrack).tracksCoordAmpCG;
    
    %get track's start time and end time
    startTime = seqOfEvents(1,1);
    endTime = seqOfEvents(end,1);
    
    %find starts of segments making this compound track
    indxStart = find( seqOfEvents(:,2)==1 & isnan(seqOfEvents(:,4)) )';
    
    %go over all starts
    for iSE = indxStart
        
        %get frame of start and track segment involved
        frameSE = seqOfEvents(iSE,1);
        segment1 = seqOfEvents(iSE,3);
        
        %make the coordinates at the start negative
        %also make the coordinate of the frame just after negative
        tracksCoordAmpCG(segment1,(frameSE-startTime)*8+1:(frameSE-startTime+2)*8) = ...
            -abs(tracksCoordAmpCG(segment1,(frameSE-startTime)*8+1:(frameSE-startTime+2)*8));
        
    end
    
    %find ends of segments making this compound track
    indxEnd = find( seqOfEvents(:,2)==2 & isnan(seqOfEvents(:,4)) )';
    
    %go over all ends
    for iSE = indxEnd
        
        %get frame of end and track segment involved
        frameSE = seqOfEvents(iSE,1);
        segment1 = seqOfEvents(iSE,3);
        
        %make the coordinates at the end negative
        %also make the coordinate of the frame just before negative
        tracksCoordAmpCG(segment1,(frameSE-startTime-1)*8+1:(frameSE-startTime+1)*8) = ...
            -abs(tracksCoordAmpCG(segment1,(frameSE-startTime-1)*8+1:(frameSE-startTime+1)*8));
        
    end
    
    %find merges and splits happening in this compound track
    indxMS = find(~isnan(seqOfEvents(:,4)))';

    %go over all merges/splits
    for iMS = indxMS

        %get frame of merge/split and track segments involved
        frameMS = seqOfEvents(iMS,1);
        segment1 = seqOfEvents(iMS,3);
        segment2 = seqOfEvents(iMS,4);
        
        numCol = size(tracksCoordAmpCG,2);

        %make the coordinates just before/after a merge/split negative
        %go 2 frames before and 2 frames after
        tracksCoordAmpCG(segment1,(frameMS-startTime-2)*8+1:min(numCol,(frameMS-startTime+2)*8)) = ...
            -abs(tracksCoordAmpCG(segment1,(frameMS-startTime-2)*8+1:min(numCol,(frameMS-startTime+2)*8)));
        tracksCoordAmpCG(segment2,(frameMS-startTime-2)*8+1:min(numCol,(frameMS-startTime+2)*8)) = ...
            -abs(tracksCoordAmpCG(segment2,(frameMS-startTime-2)*8+1:min(numCol,(frameMS-startTime+2)*8)));

    end

    %store the compound track's coordinates in big matrix
    trackedFeatureInfo(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,8*(startTime-1)+1:8*endTime) = tracksCoordAmpCG;

end

% %determine maximum coordinate values for false positive features
% xCoordMax = max(max(trackedFeatureInfo(:,1:8:end)));
% yCoordMax = max(max(trackedFeatureInfo(:,2:8:end)));

%% in each frame, delete "percentMissing+-std" of the features that can be deleted (i.e. features not just before/after a merge/split)
%% also add false positives
%% and store information in movieInfo

%initialize vector for output
percentMissingActual = zeros(numFrames,1);

%go over all frames ...
for iFrame = 1 : numFrames

    %find number of features in this frame
    numFeat = length(find(~isnan(trackedFeatureInfo(:,(iFrame-1)*8+1))));
    
    %determine number of features to delete
    numFeatDelete = round((percentMissing+randn(1)*missStd)*numFeat/100);
    
    %determine number of false positive features
    numFalsePositive = round((percentFP+randn(1)*fpStd)*numFeat/100);

    %find indices of features that can be deleted in this frame
    indxCanDelete = find(trackedFeatureInfo(:,(iFrame-1)*8+1)>0);
    
    %randomly choose which features to delete
    if length(indxCanDelete) > numFeatDelete
        deleteIndx = randsample(indxCanDelete,numFeatDelete);
    else
        deleteIndx = indxCanDelete;
    end
    
    %prevent some deletions if they generate segments between gaps that are
    %shorter than the minimum allowed segment length
    if minSegLength > 1 && ~isempty(deleteIndx)
        colIndx = iFrame-minSegLength : iFrame - 1;
        colIndx = colIndx(colIndx >= 1);
        if length(colIndx) < minSegLength
            deleteIndx = [];
        else
            tracksMatChunk = trackedFeatureInfo(deleteIndx,(colIndx-1)*8+1);
            badIndx = find(isnan(tracksMatChunk(:,1))&all(~isnan(tracksMatChunk(:,2:end)),2));
            goodIndx = setdiff(1:length(deleteIndx),badIndx);
            deleteIndx = deleteIndx(goodIndx);
        end
    end
    percentMissingActual(iFrame) = 100*length(deleteIndx)/numFeat;
    
    %delete features by making their coordinates NaN
    trackedFeatureInfo(deleteIndx,(iFrame-1)*8+1:iFrame*8) = NaN;
    
    %get the coordinates of features detected in this frame
    coord = trackedFeatureInfo(:,(iFrame-1)*8+1:(iFrame-1)*8+4);
    coord = abs(coord(~isnan(coord(:,1)),:));
    numFeat = size(coord,1);
    
    %calculate the coordinate limits for false positives
    xCoordMin = min(coord(:,1));
    xCoordMax = max(coord(:,1));
    yCoordMin = min(coord(:,2));
    yCoordMax = max(coord(:,2));
    
    %generate false positive features
    %     xCoordFP = rand(numFalsePositive,1)*xCoordMax;
    %     yCoordFP = rand(numFalsePositive,1)*yCoordMax;
    xCoordFP = rand(numFalsePositive,1)*(xCoordMax-xCoordMin)+xCoordMin;
    yCoordFP = rand(numFalsePositive,1)*(yCoordMax-yCoordMin)+yCoordMin;    
    ampFP = mean(coord(:,4))*ones(numFalsePositive,1);

    %store feature information in movieInfo
    movieInfo(iFrame).xCoord = [[coord(:,1); xCoordFP] zeros(numFeat+numFalsePositive,1)];
    movieInfo(iFrame).yCoord = [[coord(:,2); yCoordFP] zeros(numFeat+numFalsePositive,1)];
    movieInfo(iFrame).amp    = [[coord(:,4); ampFP] zeros(numFeat+numFalsePositive,1)];

end

%take average for output
percentMissingActual = mean(percentMissingActual);

%% construct tracksSimMiss which includes NaNs for "closed gaps"

%initialize tracksSimMiss
tracksSimMiss = tracksSim;

for iTrack = 1 : numTracks
    
    %get track's sequence of events
    seqOfEvents = tracksSimMiss(iTrack).seqOfEvents;

    %get track's start time and end time
    startTime = seqOfEvents(1,1);
    endTime = seqOfEvents(end,1);

    %extract track's coordinates from big matrix
    tracksCoordAmpCG = trackedFeatureInfo(trackStartRow(iTrack):...
        trackStartRow(iTrack)+numSegments(iTrack)-1,(startTime-1)*8+1:...
        endTime*8);
    tracksCoordAmpCG = abs(tracksCoordAmpCG);
    
    %store coordinates in tracksSimMiss
    tracksSimMiss(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
    
end


%% %%% ~~ the end ~~ %%%
