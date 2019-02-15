function tracks = normalizeTracks(tracks,movieInfo)
% normalizeTracks adds seqOfEvents and tracksCoordAmpCG to a tracks
% structure if they are missing
    if(~isfield(tracks,'seqOfEvents'))
        % if seqOfEvents is missing, then we assume we have a single
        % segment spanning the entire length of the movie
        [tracks.seqOfEvents] = deal([1 1 1 NaN; length(movieInfo) 2 1 NaN]);
    end
    if(~isfield(tracks,'tracksCoordAmpCG'))
        tracksCoordAmpCG = getFeatFromIdx(tracks,movieInfo);
        [tracks.tracksCoordAmpCG] = deal(tracksCoordAmpCG{:});
    end
end
