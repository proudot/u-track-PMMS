function sqDist = sqDisplacementDistr(tracks,minTrackLen)
%sqDisplacementDistr return the square displacement distribution.
%
%SYNOPSIS sqDist = sqDisplacementDistr(tracks,minTrackLen) 
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%
%OUTPUT sqDist : Distribution of square displacement.
%
%Philippe Roudot : 2013



%% input

if nargin < 1 || isempty(tracks)
    disp(' Missing input argument!');
    return
end

numTracks = length(tracks);
sqDist=[];
%go over all linear tracks ...
for iTrack = 1 : numTracks
    
    %construct matrix of linked features
    coord = tracks(iTrack).tracksCoordAmpCG;
    numSegments = size(coord,1);

    %retain only speeds of segments that are at least minTrackLen long
    segmentLen = zeros(numSegments,1);
    for iSegment = 1 : numSegments
        
        trackXVel = coord(iSegment,9:8:end) - coord(iSegment,1:8:end-8);
        trackXVel=trackXVel(coord(iSegment,9:8:end)~=0&coord(iSegment,1:8:end-8)~=0);
        trackYVel = coord(iSegment,10:8:end)- coord(iSegment,2:8:end-8);
        trackYVel=trackYVel(coord(iSegment,9:8:end)~=0&coord(iSegment,1:8:end-8)~=0);
        trackSpeed = ( trackXVel.^2 + trackYVel.^2 );
        if(size(trackXVel,2)>minTrackLen);
            sqDist = [sqDist; trackSpeed(:)];
        end 
    end
end

%% ~~~ the end ~~~
