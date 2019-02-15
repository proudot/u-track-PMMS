function tracksMat = convJCCellsToMat(tracksJC,numFrames)

numTracks = length(tracksJC);

tracksMat = NaN(numTracks,numFrames*8);

for iTrack = 1 : numTracks
    
    trackCurrent = tracksJC{iTrack};
    numRows = size(trackCurrent,1);
    
    for iRow = 1 : numRows
        
        iFrame = trackCurrent(iRow,4) + 1;
        
        iCol = (iFrame-1)*8;
        
        tracksMat(iTrack,iCol+1:iCol+3) = trackCurrent(iRow,1:3);
        
    end
    
end


