function [trackedFeatureInfo,trackedFeatureIndx] = convStruct2SparseMatNoMS(tracksFinal)
%CONVSTRUCT2SPARSEMATNOMS converts tracks from structure format to sparse matrix format, provided there are NO merges/splits.
%
%SYNPOSIS [trackedFeatureInfo,trackedFeatureIndx] = convStruct2SparseMatNoMS(tracksFinal)
%
%INPUT  tracksFinal: Output of trackCloseGapsKalman, when run with
%                    gapCloseParam.mergeSplit = 0.
%
%OUTPUT trackedFeatureInfo, trackedFeatureIndx: Output of trackWithGapClosing.
%
%REMARKS While matrix is in sparse format, gaps are still indicated by NaN
%
%Khuloud Jaqaman, August 2011

%% conversion

%get number of tracks
numTracks = length(tracksFinal);

%collect track information
trackInfo = repmat(struct('rowIndx',[],'colIndx',[],'value',[]),numTracks,1);
for iTrack = 1 : numTracks
    startTime = tracksFinal(iTrack).seqOfEvents(1,1);
    endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
    colNum = 8*(startTime-1)+1 : 8*endTime;
    trackInfo(iTrack).rowIndx = iTrack*ones(length(colNum),1);
    trackInfo(iTrack).colIndx = colNum';
    trackInfo(iTrack).value = tracksFinal(iTrack).tracksCoordAmpCG';
end

rowIndx = vertcat(trackInfo.rowIndx);
colIndx = vertcat(trackInfo.colIndx);
values = vertcat(trackInfo.value);

trackedFeatureInfo = sparse(rowIndx,colIndx,values);

if nargout == 2
    
    %collect track information
    trackInfo = repmat(struct('rowIndx',[],'colIndx',[],'value',[]),numTracks,1);
    for iTrack = 1 : numTracks
        startTime = tracksFinal(iTrack).seqOfEvents(1,1);
        endTime   = tracksFinal(iTrack).seqOfEvents(end,1);
        colNum = startTime : endTime;
        trackInfo(iTrack).rowIndx = iTrack*ones(length(colNum),1);
        trackInfo(iTrack).colIndx = colNum';
        trackInfo(iTrack).value = tracksFinal(iTrack).tracksFeatIndx';
    end
    
    rowIndx = vertcat(trackInfo.rowIndx);
    colIndx = vertcat(trackInfo.colIndx);
    values = vertcat(trackInfo.value);
    
    trackedFeatureIndx = sparse(rowIndx,colIndx,values);
    
end

%% ~~~ the end ~~~