function [fracLinksCorrect,fracLinksWrong,fracLinksFP,...
    fracGapsCorrect,fracGapsWrong,fracGapsFP,...
    fracMSCorrect,fracMSWrongAndFP,...
    aveRatioDisp2NNDist,ratioAveNNDist2MaxDisp] ...
    = summTrackPerformance(resultsSim)

n3 = size(resultsSim,3);

if n3 == 1 %simulation with no false positives
    
    [nMiss,nRun] = size(resultsSim);
    
    [fracLinksCorrect,fracLinksWrong,fracGapsCorrect,fracGapsWrong] = ...
        deal(NaN(nMiss,1));
    fracLinksFP = [];
    fracGapsFP = [];
    
    for i = 1 : nMiss
        tmp = vertcat(resultsSim(i,:).linkStats0);
        fracLinksCorrect(i) = sum(tmp(:,3))/sum(tmp(:,1));
        fracLinksWrong(i) = sum(tmp(:,4))/sum(tmp(:,1));
    end
    
    for i = 1 : nMiss
        tmp1 = vertcat(resultsSim(i,:).gapStats0);
        tmp2 = vertcat(resultsSim(i,:).gapInfo);
        fracGapsCorrect(i) = length(find(tmp1(:,2)==0))/size(tmp2,1);
        fracGapsWrong(i) = length(find(tmp1(:,2)==1))/size(tmp2,1);
    end
    
    tracks = resultsSim(1,1).tracksSimMiss;
    tracks = convStruct2MatNoMS(tracks);
    [numTrack,numFrames] = size(tracks);
    numFrames = numFrames/8;
    nnDist = NaN(numTrack,numFrames-1,nRun);
    displacement = NaN(numTrack,numFrames-1,nRun);
    for i = 1 : nRun
        tracks = resultsSim(1,i).tracksSimMiss;
        tracks = convStruct2MatNoMS(tracks);
        xCoord = tracks(:,1:8:end);
        yCoord = tracks(:,2:8:end);
        displacement(:,:,i) = sqrt((diff(xCoord,[],2)).^2 + (diff(yCoord,[],2)).^2);
        for j = 1 : numFrames-1
            indxGood = find(~isnan(xCoord(:,j)));
            nnDistFrame = createDistanceMatrix([xCoord(indxGood,j) ...
                yCoord(indxGood,j)],[xCoord(indxGood,j) yCoord(indxGood,j)]);
            nnDistFrame = sort(nnDistFrame,2);
            nnDist(indxGood,j,i) = nnDistFrame(:,2);
        end
    end
    aveRatioDisp2NNDist = displacement ./ nnDist;
    aveRatioDisp2NNDist = nanmean(aveRatioDisp2NNDist(:));
    ratioAveNNDist2MaxDisp = nanmean(nnDist(:)) / max(displacement(:));
    
else %simulation with false positives
    
    [nFP,nMiss,nRun] = size(resultsSim);
    
    [fracLinksCorrect,fracLinksWrong,fracLinksFP,...
        fracGapsCorrect,fracGapsWrong,fracGapsFP...
        fracMSCorrect,fracMSWrongAndFP] ...
        = deal(NaN(nMiss,nFP));
    
    for i = 1 : nFP
        for j = 1 : nMiss
            tmp = vertcat(resultsSim(i,j,:).linkStats0);
            fracLinksCorrect(j,i) = sum(tmp(:,3))/sum(tmp(:,1));
            fracLinksWrong(j,i) = sum(tmp(:,4))/sum(tmp(:,1));
            fracLinksFP(j,i) = sum(tmp(:,5))/sum(tmp(:,1));
        end
    end
    
    for i = 1 : nFP
        for j = 1 : nMiss
            tmp1 = vertcat(resultsSim(i,j,:).gapStats0);
            tmp2 = vertcat(resultsSim(i,j,:).gapInfo);
            fracGapsCorrect(j,i) = length(find(tmp1(:,2)==0))/size(tmp2,1);
            fracGapsWrong(j,i) = length(find(tmp1(:,2)==1))/size(tmp2,1);
            fracGapsFP(j,i) = length(find(tmp1(:,2)==2))/size(tmp2,1);
        end
    end
    
    for i = 1 : nFP
        for j = 1 : nMiss
            tmp = vertcat(resultsSim(i,j,:).mergeSplitStats0);
            fracMSCorrect(j,i) = sum(tmp(:,3))/sum(tmp(:,1));
            fracMSWrongAndFP(j,i) = (sum(tmp(:,2))-sum(tmp(:,3)))/sum(tmp(:,1));
        end
    end
    
    tracks = resultsSim(1,1,1).tracksSimMiss;
    tracks = convStruct2MatNoMS(tracks);
    [numTrack,numFrames] = size(tracks);
    numFrames = numFrames/8;
    nnDist = NaN(numTrack,numFrames-1,nRun);
    displacement = NaN(numTrack,numFrames-1,nRun);
    for i = 1 : nRun
        tracks = resultsSim(1,1,i).tracksSimMiss;
        tracks = convStruct2MatNoMS(tracks);
        xCoord = tracks(:,1:8:end);
        yCoord = tracks(:,2:8:end);
        numRow = size(xCoord,1);
        displacement(1:numRow,:,i) = sqrt((diff(xCoord,[],2)).^2 + (diff(yCoord,[],2)).^2);
        for j = 1 : numFrames-1
            indxGood = find(~isnan(xCoord(:,j)));
            nnDistFrame = createDistanceMatrix([xCoord(indxGood,j) ...
                yCoord(indxGood,j)],[xCoord(indxGood,j) yCoord(indxGood,j)]);
            nnDistFrame = sort(nnDistFrame,2);
            nnDist(indxGood,j,i) = nnDistFrame(:,2);
        end
    end
    displacement(displacement==0) = NaN;
    nnDist(nnDist==0) = NaN;
    aveRatioDisp2NNDist = displacement ./ nnDist;
    aveRatioDisp2NNDist = nanmean(aveRatioDisp2NNDist(:));
    ratioAveNNDist2MaxDisp = nanmean(nnDist(:)) / max(displacement(:));
    
end
