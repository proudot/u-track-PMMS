
%% possible links per feature

meanLinksPerFeat = NaN(8,5);
maxLinksPerFeat = NaN(8,5);
for j=1:5
    eval(['linksPerFeat = vertcat(resultsSim' num2str(j) '(:).linksPerFeat0);']);
    maxLinksPerFeatAll = max(linksPerFeat);
    for i=1:8
        eval(['linksPerFeat = vertcat(resultsSim' num2str(j) '(i,:).linksPerFeat0);']);
        eval(['freqlinksPerFeat' num2str(j) '(i,:) = hist(linksPerFeat,'...
            '(1:maxLinksPerFeatAll))/length(linksPerFeat);']);
        meanLinksPerFeat(i,j) = mean(linksPerFeat);
        maxLinksPerFeat(i,j) = max(linksPerFeat);
    end
end

clear i j linksPerFeat maxLinksPerFeatAll


%% possible links per track

meanLinksPerTrack = NaN(8,5);
maxLinksPerTrack = NaN(8,5);
for j=1:5
    eval(['linksPerTrack = vertcat(resultsSim' num2str(j) '(:).linksPerTrack0);']);
    maxLinksPerTrackAll = max(linksPerTrack);
    for i=1:8
        eval(['linksPerTrack = vertcat(resultsSim' num2str(j) '(i,:).linksPerTrack0);']);
        eval(['freqlinksPerTrack' num2str(j) '(i,:) = hist(linksPerTrack,'...
            '(1:maxLinksPerTrackAll))/length(linksPerTrack);']);
        meanLinksPerTrack(i,j) = mean(linksPerTrack);
        maxLinksPerTrack(i,j) = max(linksPerTrack);
    end
end

clear i j linksPerTrack maxLinksPerTrackAll


%% frame-to-frame displacement

%ignore merging/splitting features which are placed almost on top of their
%"other half"

meanXDisp = NaN(8,5);
meanYDisp = NaN(8,5);
meanXYDisp = NaN(8,5);
for j=1:5
    for i=1:8
        eval(['tracks = vertcat(resultsSim' num2str(j) '(i,:).tracksSimMiss);']);
        numTracks = length(tracks);
        seqOfEvents = vertcat(tracks.seqOfEvents);
        numFrames = max(seqOfEvents(:,1));
        xCoord = NaN(numTracks,numFrames);
        yCoord = xCoord;
        for m=1:numTracks
            seqOfEvents = tracks(m).seqOfEvents;
            startTime = seqOfEvents(1,1);
            endTime = seqOfEvents(end,1);
            xCoord(m,startTime:endTime) = tracks(m).tracksCoordAmpCG(1,1:8:end);
            yCoord(m,startTime:endTime) = tracks(m).tracksCoordAmpCG(1,2:8:end);
        end
        xDisp = xCoord(:,2:end) - xCoord(:,1:end-1);
        yDisp = yCoord(:,2:end) - yCoord(:,1:end-1);
        xDisp = xDisp(:);
        yDisp = yDisp(:);
        xDisp = abs(xDisp(~isnan(xDisp)));
        yDisp = abs(yDisp(~isnan(yDisp)));
        displacement(i).xDisp = xDisp;
        displacement(i).yDisp = yDisp;
        numObs(i) = length(xDisp);
        meanXDisp(i,j) = nanmean(xDisp);
        meanYDisp(i,j) = nanmean(yDisp);
        meanXYDisp(i,j) = nanmean(sqrt(xDisp.^2+yDisp.^2));
    end
    eval(['xDisplacement' num2str(j) ' = NaN(8,max(numObs));'])
    eval(['yDisplacement' num2str(j) ' = NaN(8,max(numObs));'])
    for i=1:8
        eval(['xDisplacement' num2str(j) '(i,1:numObs(i)) = displacement(i).xDisp'';'])
        eval(['yDisplacement' num2str(j) '(i,1:numObs(i)) = displacement(i).yDisp'';'])
    end
end

clear i j m tracks numTracks seqOfEvents numFrames xCoord yCoord 
clear startTime endTime xDisp yDisp displacement numObs


%% nearest-neighbor distances

%ignore merging/splitting features which are placed almost on top of their
%"other half"

meanNearestNeighborDist = NaN(8,5);
fracNnDistLessAveDisp = NaN(8,5);
fracNnDistLess2AveDisp = NaN(8,5);
for j=1:5
    for i=1:8
        nnDist = [];
        for k=1:2 %% DON'T FORGET TO CHANGE BETWEEN 5 AND 10 %%
            eval(['tracks = resultsSim' num2str(j) '(i,k).tracksSimMiss;']);
            numTracks = length(tracks);
            seqOfEvents = vertcat(tracks.seqOfEvents);
            numFrames = max(seqOfEvents(:,1));
            xCoord = NaN(numTracks,numFrames);
            yCoord = xCoord;
            for m=1:numTracks
                seqOfEvents = tracks(m).seqOfEvents;
                startTime = seqOfEvents(1,1);
                endTime = seqOfEvents(end,1);
                xCoord(m,startTime:endTime) = tracks(m).tracksCoordAmpCG(1,1:8:end);
                yCoord(m,startTime:endTime) = tracks(m).tracksCoordAmpCG(1,2:8:end);
            end
            for m=1:numFrames
                xCoord1 = xCoord(:,m);
                xCoord1 = xCoord1(~isnan(xCoord1));
                yCoord1 = yCoord(:,m);
                yCoord1 = yCoord1(~isnan(yCoord1));
                featureDist = createDistanceMatrix([xCoord1 yCoord1],...
                    [xCoord1 yCoord1]);
                featureDist = sort(featureDist,2);
                featureDist = featureDist(:,2);
                nnDist = [nnDist; featureDist];
            end
        end
        nnDist2(i).values = nnDist;
        nnDistNum(i) = length(nnDist);
        meanNearestNeighborDist(i,j) = nanmean(nnDist);
        fracNnDistLessAveDisp(i,j) = length(find(nnDist<mean(meanXYDisp(:))))/length(nnDist);
        fracNnDistLess2AveDisp(i,j) = length(find(nnDist<2*mean(meanXYDisp(:))))/length(nnDist);
    end
    eval(['nearestNeighborDist' num2str(j) ' = NaN(8,max(nnDistNum));'])
    for i=1:8
        eval(['nearestNeighborDist' num2str(j) '(i,1:nnDistNum(i)) = nnDist2(i).values'';'])
    end
end

clear i j k m numTracks seqOfEvents numFrames xCoord yCoord xCoord1 yCoord1 
clear startTime endTime featureDist nnDist nnDist2 nnDistNum tracks


