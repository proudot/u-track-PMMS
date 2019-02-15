function movieInfo = genMovieInfoFromTracks(tracksSim,percentMissing)
%GENMOVIEINFOFROMTRACKS generates a list of detected features per frame from supplied tracks
%
%SYNOPSIS movieInfo = genMovieInfoFromTracks(tracksSim,percentMissing)
%
%INPUT  tracksSim     : Output of simulateMimickCD36_MS.
%       percentMissing: Percentage of missing features in movie.
%
%OUTPUT movieInfo: List of detected features per frame, in the format
%                  required for the input of trackWithGapClosing and
%                  trackCloseGapsKalman.
%
%Khuloud Jaqaman, October 2007

%get number of frames in movie
seqOfEvents = vertcat(tracksSim.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%get number of tracks
numTracks = length(tracksSim);

%define standard deviation of missing features
missStd = 0.1*percentMissing;

%pre-allocate memory for movieInfo
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numFrames,1);

%go over all frames ...
for iFrame = 1 : numFrames
    
    %initialize variables
    xCoordNot2Delete = [];
    yCoordNot2Delete = [];
    ampNot2Delete = [];
    xCoordCanDelete = [];
    yCoordCanDelete = [];
    ampCanDelete = [];
    
    %go over all tracks ...
    for iTrack = 1 : numTracks
        
        %get track's sequence of events and coordinates
        seqOfEvents = tracksSim(iTrack).seqOfEvents;
        tracksCoordAmpCG = tracksSim(iTrack).tracksCoordAmpCG;
        
        %get track's start time and end time
        startTime = seqOfEvents(1,1);
        endTime = seqOfEvents(end,1);
        
        %if track exists in this frame ...
        if iFrame >= startTime && iFrame <= endTime

            %find indices of merges/splits happening in this frame or in
            %the next frame
            indx = find( (seqOfEvents(:,1) == iFrame | ...
                seqOfEvents(:,1) == iFrame+1) & ~isnan(seqOfEvents(:,4)) );

            %determine tracks involved
            indxTracks = [seqOfEvents(indx,3); seqOfEvents(indx,4)];
            indxTracks = unique(indxTracks);
            
            %features involved in a merge/split cannot be deleted
            xCoordNot2Delete = [xCoordNot2Delete; ...
                tracksCoordAmpCG(indxTracks,(iFrame-startTime)*8+1)];
            yCoordNot2Delete = [yCoordNot2Delete; ...
                tracksCoordAmpCG(indxTracks,(iFrame-startTime)*8+2)];
            ampNot2Delete    = [ampNot2Delete; ...
                tracksCoordAmpCG(indxTracks,(iFrame-startTime)*8+4)];
           
            %features not involved in a merge/split can be deleted
            indxTracks = setxor((1:size(tracksCoordAmpCG,1)),indxTracks);
            xCoordCanDelete = [xCoordCanDelete; ...
                tracksCoordAmpCG(indxTracks,(iFrame-startTime)*8+1)];
            yCoordCanDelete = [yCoordCanDelete; ...
                tracksCoordAmpCG(indxTracks,(iFrame-startTime)*8+2)];
            ampCanDelete    = [ampCanDelete; ...
                tracksCoordAmpCG(indxTracks,(iFrame-startTime)*8+4)];
            
        end %(if iFrame >= startTime && iFrame <= endTime)
  
    end %(for iTrack = 1 : numTracks)

    %keep only non-NaNs
    keepIndx = find(~isnan(xCoordNot2Delete));
    numFeatNot2Delete = length(keepIndx);
    xCoordNot2Delete = xCoordNot2Delete(keepIndx);
    yCoordNot2Delete = yCoordNot2Delete(keepIndx);
    ampNot2Delete    = ampNot2Delete(keepIndx);
    keepIndx = find(~isnan(xCoordCanDelete));
    numFeatCanDelete = length(keepIndx);
    xCoordCanDelete = xCoordCanDelete(keepIndx);
    yCoordCanDelete = yCoordCanDelete(keepIndx);
    ampCanDelete    = ampCanDelete(keepIndx);
    
    %get total number of features
    numFeat = numFeatNot2Delete + numFeatCanDelete;
    
    %decide on number of features to keep from those that can be deleted
    numFeatDetected = numFeatCanDelete - round( (percentMissing+randn(1)*missStd)*numFeat/100 );
    
    %choose which features to keep
    if numFeatDetected > 0
        keepIndx = randsample((1:numFeatCanDelete),numFeatDetected);
    else
        numFeatDetected = 0;
        keepIndx = [];
    end

    %store feature information in movieInfo
    movieInfo(iFrame).xCoord = [[xCoordNot2Delete; xCoordCanDelete(keepIndx)] ...
        zeros(numFeatNot2Delete+numFeatDetected,1)];
    movieInfo(iFrame).yCoord = [[yCoordNot2Delete; yCoordCanDelete(keepIndx)] ...
        zeros(numFeatNot2Delete+numFeatDetected,1)];
    movieInfo(iFrame).amp    = [[ampNot2Delete; ampCanDelete(keepIndx)] ...
        zeros(numFeatNot2Delete+numFeatDetected,1)];

end %(for iFrame = 1 : numFrames)

%%% ~~ the end ~~ %%%
