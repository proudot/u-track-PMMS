function [trackLength,gapLength,numMergesSplits] = getGapMSFractions(tracks)
%GETGAPMSFRACTIONS calculates the fraction of a track that is a gap or involves a merge or split
%
%SYNPOSIS 
%
%INPUT
%
%OUTPUT
%
%Khuloud Jaqaman, April 2008

%if track is in matrix format, convert to structure
if ~isstruct(tracks)

    %get number of tracks
    numTracks = size(tracks,1);
    
    %calculate track start times, end times and lifetimes
    trackSEL = getTrackSEL(tracks);

    %store input tracks in new variable
    tracksInput = tracks;

    %generate structure variable
    tracks = repmat(struct('tracksCoordAmpCG',[],'seqOfEvents',[]),numTracks,1);
    
    %go over all tracks and fill in the fields
    for iTrack = 1 : numTracks
        tracks(iTrack).tracksCoordAmpCG = tracksInput(iTrack,...
            8*(trackSEL(iTrack,1)-1)+1:8*trackSEL(iTrack,2));
        tracks(iTrack).seqOfEvents = [trackSEL(iTrack,1) 1 1 NaN; ...
            trackSEL(iTrack,2) 2 1 NaN];
    end

else
    
    %get number of tracks
    numTracks = length(tracks);

end

%initialize output variables
gapLength = NaN(numTracks,1);
trackLength = NaN(numTracks,1);
numMergesSplits = NaN(numTracks,1);

%go over all tracks ...
for iTrack = 1 : numTracks

    %extract the coordinates of this track
    trackCoord = tracks(iTrack).tracksCoordAmpCG;

    %get the segment start and end time information
    segmentSEL = getTrackSEL(trackCoord);
    
    %keep only the x-coordinate
    trackCoord = trackCoord(:,1:8:end);

    %go over each segment in this track
    gapLengthTmp = 0;
    trackLengthTmp = 0;
    for iSegment = 1 : size(trackCoord,1)
        
        %calculate the length of gaps in this track
        gapLengthTmp = gapLengthTmp + length(find(isnan(trackCoord(iSegment,...
            segmentSEL(iSegment,1):segmentSEL(iSegment,2)))));
        
        %calculate the length of the whole track (including gaps)
        trackLengthTmp = trackLengthTmp + segmentSEL(iSegment,3);
        
    end
    
    %store information in global array
    gapLength(iTrack) = gapLengthTmp;
    trackLength(iTrack) = trackLengthTmp;
    
    %extract the track's sequence of events
    seqOfEvents = tracks(iTrack).seqOfEvents;
    
    %calculate the number of merges and splits in track
    numMergesSplits(iTrack) = length(find(~isnan(seqOfEvents(:,4))));
    
end

%%% ~~~ the end ~~~

