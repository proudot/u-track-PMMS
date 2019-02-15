function [mergeTPFPFN,mergeEval,splitTPFPFN,splitEval] = scoreMSadvanced(...
    tracksFinal,tracksSim,srSpace,srTime,doPlot)
%SCOREMSADVANCED compares merges/splits obtained via tracking to the ground truth, returning time information as well as scores
%
%SYNOPSIS [mergeTPFPFN,mergeScore,splitTPFPFN,splitScore] = scoreMSadvanced(tracksFinal,tracksSim,srSpace,srTime)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : Simulated tracks as obtained from
%                     simulateMimickCD36_MS (structure format).
%       srSPace     : Spatial search radius to search for correspondence
%                     between merges/splits from tracker and in ground
%                     truth. Optional. Default: 5 space units.
%       srTime      : Temporal search radius to search for correspondence
%                     between merges/splits from tracker and in ground
%                     truth. Optional. Default: 5 time units.
%       doPlot      : 1 for visual output, 0 otherwise. Optional. Default: 0.
%
%OUTPUT mergeTPFPFN : Row vector storing number of true positive, false
%                     positive and false negative merges.
%       mergeEval   : Structure with detailed merge evaluation output.
%                     Contains the fields:
%           .score     : Array indicating for every merge whether it is a 
%                        true or false positive (1/0), and, if true, the 
%                        spatial and temporal distance between it and the
%                        ground truth.
%           .detailsTP : Array storing for every true positive merge the
%                        (x,y)-coordinates before and after merging in the
%                        tracking result (first layer in 3rd dimension) and
%                        in the ground truth (2nd layer). In each case,
%                        columns 1,2 are coordinates right after merging,
%                        columns 3,4 are coordinates of continued track
%                        segment right before merging, and columns 5,6 are
%                        coordinates of terminated track segment right
%                        before merging.
%           .detailsFP : Same as above, but for false positive merges.
%                        No 3rd dimension; by definition there
%                        is no ground truth.
%           .detailsFN : Same as above, but for false negative merges. 
%                        No 3rd dimension; by definition there
%                        are no tracking results.
%       splitTPFPFN : Row vector storing number of true positive, false
%                     positive and false negative splits.
%       splitEval   : Structure with detailed split evaluation output.
%                     Equivalent to mergeEval.
%
%Khuloud Jaqaman, January 2015

%% input

if nargin < 2
    error('scoreMSadvanced: Please supply tracks and ground truth');
end

if nargin < 3 || isempty(srSpace)
    srSpace = 5;
end

if nargin < 4 || isempty(srTime)
    srTime = 5;
end

if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end

%% process input variables

%get merges and splits in tracking results and ground truth
[mergesInfo0,splitsInfo0,mergesInfoSpace0,splitsInfoSpace0,...
    mergesInfoSpaceExt0,splitsInfoSpaceExt0,mergesInfoAmp0,splitsInfoAmp0] ...
    = findMergesSplits(tracksSim,2,0,0,0);
[mergesInfoF,splitsInfoF,mergesInfoSpaceF,splitsInfoSpaceF,...
    mergesInfoSpaceExtF,splitsInfoSpaceExtF,mergesInfoAmpF,splitsInfoAmpF] ...
    = findMergesSplits(tracksFinal,2,0,0,0);

%reformat merge and split information
if ~isempty(mergesInfo0)
    numMerges0 = sum(mergesInfo0(:,3));
    merges0 = [reshape(mergesInfo0(:,4:end),[],1) reshape(mergesInfoSpace0(:,1:2:end),[],1) reshape(mergesInfoSpace0(:,2:2:end),[],1)];
    mergeDetails0 = [merges0(:,2:3) reshape(mergesInfoSpaceExt0(:,1:2:end,1),[],1) reshape(mergesInfoSpaceExt0(:,2:2:end,1),[],1) ...
        reshape(mergesInfoSpaceExt0(:,1:2:end,2),[],1) reshape(mergesInfoSpaceExt0(:,2:2:end,2),[],1)];
    mergeDetailsAmp0 = [reshape(mergesInfoAmp0(:,1:3:end),[],1) reshape(mergesInfoAmp0(:,2:3:end),[],1) reshape(mergesInfoAmp0(:,3:3:end),[],1)];
    mergeTrackIndx0 = repmat(mergesInfo0(:,1),max(mergesInfo0(:,3)),1);
    indxKeep = find(merges0(:,1)~=0);
    merges0 = merges0(indxKeep,:);
    mergeDetails0 = mergeDetails0(indxKeep,:);
    mergeDetailsAmp0 = mergeDetailsAmp0(indxKeep,:);
    mergeTrackIndx0 = mergeTrackIndx0(indxKeep);
else
    numMerges0 = 0;
end
if ~isempty(splitsInfo0)
    numSplits0 = sum(splitsInfo0(:,3));
    splits0 = [reshape(splitsInfo0(:,4:end),[],1) reshape(splitsInfoSpace0(:,1:2:end),[],1) reshape(splitsInfoSpace0(:,2:2:end),[],1)];
    splitDetails0 = [splits0(:,2:3) reshape(splitsInfoSpaceExt0(:,1:2:end,1),[],1) reshape(splitsInfoSpaceExt0(:,2:2:end,1),[],1) ...
        reshape(splitsInfoSpaceExt0(:,1:2:end,2),[],1) reshape(splitsInfoSpaceExt0(:,2:2:end,2),[],1)];
    splitDetailsAmp0 = [reshape(splitsInfoAmp0(:,1:3:end),[],1) reshape(splitsInfoAmp0(:,2:3:end),[],1) reshape(splitsInfoAmp0(:,3:3:end),[],1)];
    splitTrackIndx0 = repmat(splitsInfo0(:,1),max(splitsInfo0(:,3)),1);
    indxKeep = find(splits0(:,1)~=0);
    splits0 = splits0(indxKeep,:);
    splitDetails0 = splitDetails0(indxKeep,:);
    splitDetailsAmp0 = splitDetailsAmp0(indxKeep,:);
    splitTrackIndx0 = splitTrackIndx0(indxKeep);
else
    numSplits0 = 0;
end
if ~isempty(mergesInfoF)
    numMergesF = sum(mergesInfoF(:,3));
    mergesF = [reshape(mergesInfoF(:,4:end),[],1) reshape(mergesInfoSpaceF(:,1:2:end),[],1) reshape(mergesInfoSpaceF(:,2:2:end),[],1)];
    mergeDetailsF = [mergesF(:,2:3) reshape(mergesInfoSpaceExtF(:,1:2:end,1),[],1) reshape(mergesInfoSpaceExtF(:,2:2:end,1),[],1) ...
        reshape(mergesInfoSpaceExtF(:,1:2:end,2),[],1) reshape(mergesInfoSpaceExtF(:,2:2:end,2),[],1)];
    mergeDetailsAmpF = [reshape(mergesInfoAmpF(:,1:3:end),[],1) reshape(mergesInfoAmpF(:,2:3:end),[],1) reshape(mergesInfoAmpF(:,3:3:end),[],1)];
    mergeTrackIndxF = repmat(mergesInfoF(:,1),max(mergesInfoF(:,3)),1);
    indxKeep = find(mergesF(:,1)~=0);
    mergesF = mergesF(indxKeep,:);
    mergeDetailsF = mergeDetailsF(indxKeep,:);
    mergeDetailsAmpF = mergeDetailsAmpF(indxKeep,:);
    mergeTrackIndxF = mergeTrackIndxF(indxKeep);
else
    numMergesF = 0;
end
if ~isempty(splitsInfoF)
    numSplitsF = sum(splitsInfoF(:,3));
    splitsF = [reshape(splitsInfoF(:,4:end),[],1) reshape(splitsInfoSpaceF(:,1:2:end),[],1) reshape(splitsInfoSpaceF(:,2:2:end),[],1)];
    splitDetailsF = [splitsF(:,2:3) reshape(splitsInfoSpaceExtF(:,1:2:end,1),[],1) reshape(splitsInfoSpaceExtF(:,2:2:end,1),[],1) ...
        reshape(splitsInfoSpaceExtF(:,1:2:end,2),[],1) reshape(splitsInfoSpaceExtF(:,2:2:end,2),[],1)];
    splitDetailsAmpF = [reshape(splitsInfoAmpF(:,1:3:end),[],1) reshape(splitsInfoAmpF(:,2:3:end),[],1) reshape(splitsInfoAmpF(:,3:3:end),[],1)];
    splitTrackIndxF = repmat(splitsInfoF(:,1),max(splitsInfoF(:,3)),1);
    indxKeep = find(splitsF(:,1)~=0);
    splitsF = splitsF(indxKeep,:);
    splitDetailsF = splitDetailsF(indxKeep,:);
    splitDetailsAmpF = splitDetailsAmpF(indxKeep,:);
    splitTrackIndxF = splitTrackIndxF(indxKeep);
else
    numSplitsF = 0;
end

%% score merges

if numMerges0 > 0
    
    if numMergesF > 0
        
        %calculate spatial and temporal distances between tracking result and
        %ground truth
        spaceMat = createDistanceMatrix(mergesF(:,2:3),merges0(:,2:3));
        timeMat = -createDistanceMatrix(mergesF(:,1),merges0(:,1));
        
        %impose upper limits
        spaceMat(spaceMat > srSpace) = NaN;
        timeMat(abs(timeMat) > srTime) = NaN;
        
        %calculate linking cost matrix
        costMat = spaceMat + abs(timeMat);
        costMat(isnan(costMat)) = -1;
        
        %match via LAP
        linkFto0 = double(lap(costMat,-1,0,1));
        linkFto0 = linkFto0(1:numMergesF);
        
        %evaluate merges
        indxGood = find(linkFto0 <= numMerges0);
        indxFP = setdiff(1:numMergesF,indxGood);
        indxMatch = linkFto0(indxGood);
        indxFN = setdiff((1:numMerges0)',indxMatch);
        indxLin = sub2ind([numMergesF numMerges0],indxGood,indxMatch);
        mergeTP = length(indxGood);
        mergeFP = numMergesF - mergeTP;
        mergeFN = numMerges0 - mergeTP;
        mergeTPFPFN = [mergeTP mergeFP mergeFN];
        
        %detailed information for output
        mergeScore = [timeMat(indxLin) spaceMat(indxLin)];
        mergeDetailsTP(:,:,1) = mergeDetailsF(indxGood,:);
        mergeDetailsTP(:,:,2) = mergeDetails0(indxMatch,:);
        mergeDetailsFP = mergeDetailsF(indxFP,:);
        mergeDetailsFN = mergeDetails0(indxFN,:);
        mergeDetailsAmpTP(:,:,1) = mergeDetailsAmpF(indxGood,:);
        mergeDetailsAmpTP(:,:,2) = mergeDetailsAmp0(indxMatch,:);
        mergeDetailsAmpFP = mergeDetailsAmpF(indxFP,:);
        mergeDetailsAmpFN = mergeDetailsAmp0(indxFN,:);
        mergeTimeTP = [mergesF(indxGood,1) merges0(indxMatch,1)];
        mergeTimeFP = mergesF(indxFP,1);
        mergeTimeFN = merges0(indxFN,1);
        mergeTrackIndxTP = [mergeTrackIndxF(indxGood) mergeTrackIndx0(indxMatch)];
        mergeTrackIndxFP = mergeTrackIndxF(indxFP);
        mergeTrackIndxFN = mergeTrackIndx0(indxFN);
        
    else
        
        mergeScore = [];
        mergeTPFPFN = [0 0 numMerges0];
        mergeDetailsTP = [];
        mergeDetailsFP = [];
        mergeDetailsFN = mergeDetails0;
        mergeDetailsAmpTP = [];
        mergeDetailsAmpFP = [];
        mergeDetailsAmpFN = mergeDetailsAmp0;
        mergeTimeTP = [];
        mergeTimeFP = [];
        mergeTimeFN = merges0(:,1);
        mergeTrackIndxTP = [];
        mergeTrackIndxFP = [];
        mergeTrackIndxFN = mergeTrackIndx0;
        
    end
    
else
    
    mergeScore = [];
    mergeTPFPFN = [0 numMergesF 0];
    mergeDetailsTP = [];
    mergeDetailsFP = mergeDetailsF;
    mergeDetailsFN = [];
    mergeDetailsAmpTP = [];
    mergeDetailsAmpFP = mergeDetailsAmpF;
    mergeDetailsAmpFN = [];
    mergeTimeTP = [];
    mergeTimeFP = mergesF(:,1);
    mergeTimeFN = [];
    mergeTrackIndxTP = [];
    mergeTrackIndxFP = mergeTrackIndxF;
    mergeTrackIndxFN = [];
    
end

detailsTP = struct('trackIndx',mergeTrackIndxTP,'time',mergeTimeTP,'coord',mergeDetailsTP,'amp',mergeDetailsAmpTP,'score',mergeScore);
detailsFP = struct('trackIndx',mergeTrackIndxFP,'time',mergeTimeFP,'coord',mergeDetailsFP,'amp',mergeDetailsAmpFP);
detailsFN = struct('trackIndx',mergeTrackIndxFN,'time',mergeTimeFN,'coord',mergeDetailsFN,'amp',mergeDetailsAmpFN);
mergeEval = struct('detailsTP',detailsTP,'detailsFP',detailsFP,'detailsFN',detailsFN);

%% score splits

if numSplits0 > 0
    
    if numSplitsF > 0
        
        %calculate spatial and temporal distances between tracking result and
        %ground truth
        spaceMat = createDistanceMatrix(splitsF(:,2:3),splits0(:,2:3));
        timeMat = -createDistanceMatrix(splitsF(:,1),splits0(:,1));
        
        %impose upper limits
        spaceMat(spaceMat > 5) = NaN;
        timeMat(abs(timeMat) > 5) = NaN;
        
        %calculate linking cost matrix
        costMat = spaceMat + abs(timeMat);
        costMat(isnan(costMat)) = -1;
        
        %match via LAP
        linkFto0 = double(lap(costMat,-1,0,1));
        linkFto0 = linkFto0(1:numSplitsF);
        
        %evaluate splits
        indxGood = find(linkFto0 <= numSplits0);
        indxFP = setdiff(1:numSplitsF,indxGood);
        indxMatch = linkFto0(indxGood);
        indxFN = setdiff((1:numSplits0)',indxMatch);
        indxLin = sub2ind([numSplitsF numSplits0],indxGood,indxMatch);
        splitTP = length(indxGood);
        splitFP = numSplitsF - splitTP;
        splitFN = numSplits0 - splitTP;
        splitTPFPFN = [splitTP splitFP splitFN];
        
        %detailed information for output
        splitScore = [timeMat(indxLin) spaceMat(indxLin)];
        splitDetailsTP(:,:,1) = splitDetailsF(indxGood,:);
        splitDetailsTP(:,:,2) = splitDetails0(indxMatch,:);
        splitDetailsFP = splitDetailsF(indxFP,:);
        splitDetailsFN = splitDetails0(indxFN,:);
        splitDetailsAmpTP(:,:,1) = splitDetailsAmpF(indxGood,:);
        splitDetailsAmpTP(:,:,2) = splitDetailsAmp0(indxMatch,:);
        splitDetailsAmpFP = splitDetailsAmpF(indxFP,:);
        splitDetailsAmpFN = splitDetailsAmp0(indxFN,:);
        splitTimeTP = [splitsF(indxGood,1) splits0(indxMatch,1)];
        splitTimeFP = splitsF(indxFP,1);
        splitTimeFN = splits0(indxFN,1);
        splitTrackIndxTP = [splitTrackIndxF(indxGood) splitTrackIndx0(indxMatch)];
        splitTrackIndxFP = splitTrackIndxF(indxFP);
        splitTrackIndxFN = splitTrackIndx0(indxFN);
        
    else
        
        splitScore = [];
        splitTPFPFN = [0 0 numSplits0];
        splitDetailsTP = [];
        splitDetailsFP = [];
        splitDetailsFN = splitDetails0;
        splitDetailsAmpTP = [];
        splitDetailsAmpFP = [];
        splitDetailsAmpFN = splitDetails0;
        splitTimeTP = [];
        splitTimeFP = [];
        splitTimeFN = splits0(:,1);
        splitTrackIndxTP = [];
        splitTrackIndxFP = [];
        splitTrackIndxFN = splitTrackIndx0;
        
    end
    
else
    
    splitScore = [];
    splitTPFPFN = [0 numSplitsF 0];
    splitDetailsTP = [];
    splitDetailsFP = splitDetailsF;
    splitDetailsFN = [];
    splitDetailsAmpTP = [];
    splitDetailsAmpFP = splitDetailsF;
    splitDetailsAmpFN = [];
    splitTimeTP = [];
    splitTimeFP = splitsF(:,1);
    splitTimeFN = [];
    splitTrackIndxTP = [];
    splitTrackIndxFP = splitTrackIndxF;
    splitTrackIndxFN = [];
 
end

detailsTP = struct('trackIndx',splitTrackIndxTP,'time',splitTimeTP,'coord',splitDetailsTP,'amp',splitDetailsAmpTP,'score',splitScore);
detailsFP = struct('trackIndx',splitTrackIndxFP,'time',splitTimeFP,'coord',splitDetailsFP,'amp',splitDetailsAmpFP);
detailsFN = struct('trackIndx',splitTrackIndxFN,'time',splitTimeFN,'coord',splitDetailsFN,'amp',splitDetailsAmpFN);
splitEval = struct('detailsTP',detailsTP,'detailsFP',detailsFP,'detailsFN',detailsFN);

%% plot

if doPlot
    
    maxCoord = ceil(max([mergeDetails0(:); mergeDetailsF(:); splitDetails0(:); splitDetailsF(:)]))+1;
    
    %merges
    figure('Name','Merges')
    imshow(ones(maxCoord));
    axH = gca;
    set(axH,'visible','on');
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');
    if ~isempty(mergeDetailsTP)
        lineWithGaps(mergeDetailsTP(:,[1 3],2)',mergeDetailsTP(:,[2 4],2)',[],'Color','k');
        lineWithGaps(mergeDetailsTP(:,[1 3],1)',mergeDetailsTP(:,[2 4],1)',[],'Color','g');
        lineWithGaps(squeeze(mergeDetailsTP(:,1,:))',squeeze(mergeDetailsTP(:,2,:))',[],'Color',[0.7 0.7 0.7],'LineStyle','--');
    end
    if ~isempty(mergeDetailsFP)
        lineWithGaps(mergeDetailsFP(:,[1 3])',mergeDetailsFP(:,[2 4])',[],'Color','r');
    end
    if ~isempty(mergeDetailsFN)
        lineWithGaps(mergeDetailsFN(:,[1 3])',mergeDetailsFN(:,[2 4])',[],'Color','b');
    end
    legend('TP GT','TP','Match','FP','FN')
    if ~isempty(mergeDetailsTP)
        lineWithGaps(mergeDetailsTP(:,[1 3],2)',mergeDetailsTP(:,[2 4],2)',[],'Color','k','Marker','.');
        lineWithGaps(mergeDetailsTP(:,[1 5],2)',mergeDetailsTP(:,[2 6],2)',[],'Color','k','Marker','.','LineStyle',':');
        lineWithGaps(mergeDetailsTP(:,[1 3],1)',mergeDetailsTP(:,[2 4],1)',[],'Color','g','Marker','.');
        lineWithGaps(mergeDetailsTP(:,[1 5],1)',mergeDetailsTP(:,[2 6],1)',[],'Color','g','Marker','.','LineStyle',':');
        lineWithGaps(squeeze(mergeDetailsTP(:,1,:))',squeeze(mergeDetailsTP(:,2,:))',[],'Color',[0.7 0.7 0.7],'LineStyle','--');
    end
    if ~isempty(mergeDetailsFP)
        lineWithGaps(mergeDetailsFP(:,[1 3])',mergeDetailsFP(:,[2 4])',[],'Color','r','Marker','.');
        lineWithGaps(mergeDetailsFP(:,[1 5])',mergeDetailsFP(:,[2 6])',[],'Color','r','Marker','.','LineStyle',':');
    end
    if ~isempty(mergeDetailsFN)
        lineWithGaps(mergeDetailsFN(:,[1 3])',mergeDetailsFN(:,[2 4])',[],'Color','b','Marker','.');
        lineWithGaps(mergeDetailsFN(:,[1 5])',mergeDetailsFN(:,[2 6])',[],'Color','b','Marker','.','LineStyle',':');
    end
    
    %splits
    figure('Name','Splits')
    imshow(ones(maxCoord));
    axH = gca;
    set(axH,'visible','on');
    xlabel('x-coordinate (pixels)');
    ylabel('y-coordinate (pixels)');
    if ~isempty(splitDetailsTP)
        lineWithGaps(splitDetailsTP(:,[1 3],2)',splitDetailsTP(:,[2 4],2)',[],'Color','k');
        lineWithGaps(splitDetailsTP(:,[1 3],1)',splitDetailsTP(:,[2 4],1)',[],'Color','g');
        lineWithGaps(squeeze(splitDetailsTP(:,1,:))',squeeze(splitDetailsTP(:,2,:))',[],'Color',[0.7 0.7 0.7],'LineStyle','--');
    end
    if ~isempty(splitDetailsFP)
        lineWithGaps(splitDetailsFP(:,[1 3])',splitDetailsFP(:,[2 4])',[],'Color','r');
    end
    if ~isempty(splitDetailsFN)
        lineWithGaps(splitDetailsFN(:,[1 3])',splitDetailsFN(:,[2 4])',[],'Color','b');
    end
    legend('TP GT','TP','Match','FP','FN')
    if ~isempty(splitDetailsTP)
        lineWithGaps(splitDetailsTP(:,[1 3],2)',splitDetailsTP(:,[2 4],2)',[],'Color','k','Marker','.');
        lineWithGaps(splitDetailsTP(:,[1 5],2)',splitDetailsTP(:,[2 6],2)',[],'Color','k','Marker','.','LineStyle',':');
        lineWithGaps(splitDetailsTP(:,[1 3],1)',splitDetailsTP(:,[2 4],1)',[],'Color','g','Marker','.');
        lineWithGaps(splitDetailsTP(:,[1 5],1)',splitDetailsTP(:,[2 6],1)',[],'Color','g','Marker','.','LineStyle',':');
        lineWithGaps(squeeze(splitDetailsTP(:,1,:))',squeeze(splitDetailsTP(:,2,:))',[],'Color',[0.7 0.7 0.7],'LineStyle','--');
    end
    if ~isempty(splitDetailsFP)
        lineWithGaps(splitDetailsFP(:,[1 3])',splitDetailsFP(:,[2 4])',[],'Color','r','Marker','.');
        lineWithGaps(splitDetailsFP(:,[1 5])',splitDetailsFP(:,[2 6])',[],'Color','r','Marker','.','LineStyle',':');
    end
    if ~isempty(splitDetailsFN)
        lineWithGaps(splitDetailsFN(:,[1 3])',splitDetailsFN(:,[2 4])',[],'Color','b','Marker','.');
        lineWithGaps(splitDetailsFN(:,[1 5])',splitDetailsFN(:,[2 6])',[],'Color','b','Marker','.','LineStyle',':');
    end
    
end

