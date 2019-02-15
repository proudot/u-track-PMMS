function tracksFinal = convertMat2Struct2(tracksMat)

%get number of tracks
numTracks = size(tracksMat,1);

%reserve memory for structure storing tracks
tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
    'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracks,1);

%get the start and end time of tracks
trackSEL = getTrackSEL(tracksMat);

%go over all tracks and store information
for iTrack = 1 : numTracks
   
    %track start time and end time
    startTime = trackSEL(iTrack,1);
    endTime = trackSEL(iTrack,2);
    
    %feature indices
    %     tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxLink(iTrack,startTime:endTime);
    
    %feature coordinates and amplitudes
    tracksFinal(iTrack).tracksCoordAmpCG = full(tracksMat(iTrack,...
        (startTime-1)*8+1:endTime*8));
    
    %sequence of events
    tracksFinal(iTrack).seqOfEvents = [startTime 1 1 NaN; endTime 2 1 NaN];
    
end
