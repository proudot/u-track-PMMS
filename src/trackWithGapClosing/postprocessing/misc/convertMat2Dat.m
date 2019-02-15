function convertMat2Dat(tracksFinal,directory2save,interpolateGap)
%CONVERTMAT2DAT write the output of trackCloseGapsKalman into text files for the Kusumi lab
%
%SYNPOSIS convertMat2Dat(tracksFinal,directory2save,interpolateGap)
%
%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       directory2save: Directory where text files are to be saved.
%       interpolateGap: 1 to interpolate coordinates for closed gaps, 0 to
%                       use value from last detection
%
%REMARKS Only tracks without merging and splitting, tracks exhibiting 
%one of the following scenarios: (1) single merge, (2) single split, 
%(3) split followed by merge, and (4) merge followed by split, are
%converted
%
%Khuloud Jaqaman, February 2009

%resolve simple merges and splits
tracksFinal = resolveMergesSplitsSimple(tracksFinal);

%get number of tracks
numTracks = length(tracksFinal);

%get number of segments in each track
numSegments = ones(numTracks,1);
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
end

%keep only tracks with one segment
tracksFinal = tracksFinal(numSegments==1);

%get number of surviving tracks 
numTracks = length(tracksFinal);

%go over the tracks and store each in its file
for iTrack = 1 : numTracks
    
    %get current track information
    currentTrack = tracksFinal(iTrack);
    frameNum = (currentTrack.seqOfEvents(1,1):currentTrack.seqOfEvents(2,1))';
    xCoord = currentTrack.tracksCoordAmpCG(1:8:end)';
    yCoord = currentTrack.tracksCoordAmpCG(2:8:end)';
    amp = currentTrack.tracksCoordAmpCG(4:8:end)';
    ampDiff = [0; diff(amp)];
    
    %give frames a flag of 0 if there is a detection, and a flag of 1 if
    %there is no detection
    frameFlag = isnan(xCoord);
    goodFrames = find(~frameFlag);
    gapFrames = find(frameFlag);
    numGapFrames = length(gapFrames);
    
    %fill in missing information in gap frames, based on input flag
    if interpolateGap %interpolation

        %for every gap frame, find good frame just before it and good frame
        %just after it
        goodFrameBefore = ones(numGapFrames,1);
        goodFrameAfter = ones(numGapFrames,1);
        for iGapFrame = 1 : numGapFrames
            goodFrameBefore(iGapFrame) = goodFrames(find(goodFrames<gapFrames(iGapFrame),1,'last'));
            goodFrameAfter(iGapFrame) = goodFrames(find(goodFrames>gapFrames(iGapFrame),1,'first'));
        end

        %interpolate the value between the two good frames before and after
        xCoordChange = xCoord(goodFrameAfter) - xCoord(goodFrameBefore);
        yCoordChange = yCoord(goodFrameAfter) - yCoord(goodFrameBefore);
        frameNumIncrement = gapFrames - goodFrameBefore;
        frameNumChange = goodFrameAfter - goodFrameBefore;
        frameNumRelInc = frameNumIncrement ./ frameNumChange;
        xCoord(gapFrames) = xCoordChange .* frameNumRelInc + xCoord(goodFrameBefore);
        yCoord(gapFrames) = yCoordChange .* frameNumRelInc + yCoord(goodFrameBefore);
        
    else %no interpolation

        %for every gap frame, find good frame just before it
        correspondingGoodFrame = ones(numGapFrames,1);
        for iGapFrame = 1 : numGapFrames
            correspondingGoodFrame(iGapFrame) = goodFrames(find(goodFrames<gapFrames(iGapFrame),1,'last'));
        end

        %use the value in that frame to fill closed gaps
        xCoord(gapFrames) = xCoord(correspondingGoodFrame);
        yCoord(gapFrames) = yCoord(correspondingGoodFrame);
        
    end
        
    amp(isnan(amp)) = 0;
    ampDiff(isnan(ampDiff)) = 0;
    
    %open the file to save in
    filename2save = ['track' num2str(iTrack) '.dat'];
    fid = fopen(fullfile(directory2save,filename2save),'w+');
    
    %write information in text file
    for iFrame = 1 : length(frameNum)
        if frameFlag(iFrame) == 0
            fprintf(fid,'%4d %10.4f %10.4f OK  %7.4f %7.4f \n',frameNum(iFrame),xCoord(iFrame),yCoord(iFrame),amp(iFrame),ampDiff(iFrame));
        else
            fprintf(fid,'%4d %10.4f %10.4f Cmp %7.4f %7.4f \n',frameNum(iFrame),xCoord(iFrame),yCoord(iFrame),amp(iFrame),ampDiff(iFrame));
        end
    end
    
    %close file
    fclose(fid);
    
end
