function tracksFinal = convMetamorphTracks2Struct(metamorphTracks)

%get number of tracks
numTracks = max(metamorphTracks(:,2));

%get number of frames
numFrames = max(metamorphTracks(:,1));

%reserve memory for x and y coordinates
xCoordMat = NaN(numTracks,numFrames);
yCoordMat = xCoordMat;

%go over all tracks
for iTrack = 1 : numTracks
   
    %get the current tracks' information
    trackInfo = metamorphTracks(metamorphTracks(:,2)==iTrack,:);
    
    %get frames where current track exists
    framesExist = trackInfo(:,1);
    
    %store the track's coordinates
    xCoordMat(iTrack,framesExist) = trackInfo(:,3);
    yCoordMat(iTrack,framesExist) = trackInfo(:,4);
    
end

%assign an amplitude of 1 wherever there are coordinates
ampMat = double(~isnan(xCoordMat));
ampMat(ampMat == 0) = NaN;

%put coordinates into a big matrix of the same style as the output of
%trackWithGapClosing
tracksMat = zeros(numTracks,numFrames*8);
tracksMat(:,1:8:end) = xCoordMat;
tracksMat(:,2:8:end) = yCoordMat;
tracksMat(:,4:8:end) = ampMat;

%get the track start and end times
trackSEL = getTrackSEL(tracksMat);

%convert track matrix to track structure, i.e. the output of
%trackCloseGapsKalman
tracksFinal = repmat(struct('tracksCoordAmpCG',[],'seqOfEvents',[]),numTracks,1);
for iTrack = 1 : numTracks
    startTime = trackSEL(iTrack,1);
    endTime = trackSEL(iTrack,2);
    tracksFinal(iTrack).tracksCoordAmpCG = tracksMat(iTrack,(startTime-1)*8+1:endTime*8);
    tracksFinal(iTrack).seqOfEvents = [startTime 1 1 NaN; endTime 2 1 NaN];
end