function [analysisResPart,analysisResFull] = getMotionMergeSplitStatsAllIndBS(...
    trackData,diffAnalysisParam,msStatsParam,bsSampleSize)

%% Output
analysisResPart = [];
analysisResFull = [];

%% Input

if nargin < 3
    disp('Input arguments missing');
    return
end
   
%diffusion analysis parameters
extractType = diffAnalysisParam.extractType;
probDim = diffAnalysisParam.probDim;
checkAsym = diffAnalysisParam.checkAsym;
alphaValues = diffAnalysisParam.alphaValues;
confRadMin = diffAnalysisParam.confRadMin;

%merge/split analysis parameters
minTrackLen = msStatsParam.minTrackLen;
removePotArtifacts = msStatsParam.removePotArtifacts;

%get number of cells
numCells = length(trackData);

%put all tracks together
tracksAll = vertcat(trackData.tracks);

%% Perform overall diffusion and merge/split analyses

%diffusion analysis
diffAnalysisAll = trackDiffusionAnalysis1(tracksAll,extractType,probDim,...
    checkAsym,alphaValues,0,confRadMin);

%diffusion analysis summary
[motionStatsAll.probMotionType,motionStatsAll.motionChar] = ...
    summarizeDiffAnRes(tracksAll,minTrackLen,probDim,diffAnalysisAll,extractType);

%merge/split analysis
[msStatsAll.statsGeneral,msStatsAll.motionTypeFrac,...
    msStatsAll.probMSIfTypeFull,msStatsAll.probMSIfType12,...
    msStatsAll.numEventsPerCatFull,msStatsAll.motionRedTypeFrac,...
    msStatsAll.probMSIfRedTypeFull,msStatsAll.probMSIfRedType12] = ...
    calcStatsMS(tracksAll,minTrackLen,probDim,diffAnalysisAll,...
    removePotArtifacts);


%% Generate overall bootstrap sample

%initialize variables where results will be stored
motionStatsAllSample = repmat(struct('probMotionType',[],'motionChar',[]),...
    bsSampleSize,1);
msStatsAllSample = repmat(struct('statsGeneral',[],'motionTypeFrac',[],...
    'probMSIfTypeFull',[],'probMSIfType12',[],'numEventsPerCatFull',[],...
    'motionRedTypeFrac',[],'probMSIfRedTypeFull',[],'probMSIfRedType12',[]),...
    bsSampleSize,1);

%put the original overall statistics as the first entry
motionStatsAllSample(1) = motionStatsAll;
msStatsAllSample(1) = msStatsAll;

%remove tracks shorter than minTrackLen
criteria.lifeTime.min = minTrackLen;
indxKeep = chooseTracks(tracksAll,criteria);
tracksTmp = tracksAll(indxKeep);
diffAnalysisTmp = diffAnalysisAll(indxKeep);

%get number of remaining overall tracks
numTracks = length(tracksTmp);

%generate (bsSampleSize - 1) bootstrap samples
for iBoot = 2 : bsSampleSize
    
    %select a random subset of tracks with replacement
    randIndx = randsample((1:numTracks),numTracks,'true');
    tracksSample = tracksTmp(randIndx);
    diffAnalysisSample = diffAnalysisTmp(randIndx);
    
    %get motion statistics
    [motionStatsAllSample(iBoot).probMotionType,motionStatsAllSample(iBoot).motionChar] = ...
        summarizeDiffAnRes(tracksSample,minTrackLen,probDim,...
        diffAnalysisSample,extractType);
    
    %get merge/split statistics
    [msStatsAllSample(iBoot).statsGeneral,msStatsAllSample(iBoot).motionTypeFrac,...
        msStatsAllSample(iBoot).probMSIfTypeFull,msStatsAllSample(iBoot).probMSIfType12,...
        msStatsAllSample(iBoot).numEventsPerCatFull,msStatsAllSample(iBoot).motionRedTypeFrac,...
        msStatsAllSample(iBoot).probMSIfRedTypeFull,msStatsAllSample(iBoot).probMSIfRedType12] = ...
        calcStatsMS(tracksSample,minTrackLen,probDim,diffAnalysisSample,removePotArtifacts);
    
end

%% Perform single cell diffusion and merge/split analyses

%go over all cells
for iCell = numCells : -1 : 1
    
    %diffusion analysis
    diffAnalysisInd(iCell).results = trackDiffusionAnalysis1(...
        trackData(iCell).tracks,extractType,probDim,checkAsym,alphaValues,...
        0,confRadMin);
    
    %diffusion analysis summary
    [motionStatsInd(iCell).probMotionType,motionStatsInd(iCell).motionChar] = ...
        summarizeDiffAnRes(trackData(iCell).tracks,minTrackLen,probDim,...
        diffAnalysisInd(iCell).results,extractType);
    
    %merge/split analysis
    [msStatsInd(iCell).statsGeneral,msStatsInd(iCell).motionTypeFrac,...
        msStatsInd(iCell).probMSIfTypeFull,msStatsInd(iCell).probMSIfType12,...
        msStatsInd(iCell).numEventsPerCatFull,msStatsInd(iCell).motionRedTypeFrac,...
        msStatsInd(iCell).probMSIfRedTypeFull,msStatsInd(iCell).probMSIfRedType12] = ...
        calcStatsMS(trackData(iCell).tracks,minTrackLen,probDim,...
        diffAnalysisInd(iCell).results,removePotArtifacts);
end

%% Calculate number of repetitions for each cell to make bootstrap sample

%get number of features per cell
numFeatPerFrameInd = vertcat(msStatsInd.statsGeneral);
numFeatPerFrameInd = numFeatPerFrameInd(:,1);

%square number of features per cell
numFeatPerFrameIndSquared = numFeatPerFrameInd.^2;

%calculate relative weight for motion statistics
relWeightIndFromNumFeat = numFeatPerFrameInd/min(numFeatPerFrameInd);

%calculate relative weight for merge/split statistics
relWeightIndFromNumFeatSquared = numFeatPerFrameIndSquared/min(numFeatPerFrameIndSquared);

%get sample size resulting from features and squared features
sampSizeFeat = sum(relWeightIndFromNumFeat);
sampSizeSquaredFeat = sum(relWeightIndFromNumFeatSquared);

%divide minimum bootstrap sample size by actual sample size to get
%necessary number of repetitions to achieve minimum sample size
numRepFeat = max([1 bsSampleSize/sampSizeFeat]);
numRepSquaredFeat = max([1 bsSampleSize/sampSizeSquaredFeat]);

%calculate number of repetitions of each cell
numRepIndFromNumFeat = round(relWeightIndFromNumFeat*numRepFeat);
numRepIndFromNumFeatSquared = round(relWeightIndFromNumFeatSquared*numRepSquaredFeat);

%% Generate bootstrap samples per cell based on the number of repetitions

%determine largest number of repetitions
maxNumRep = max([numRepIndFromNumFeatSquared; numRepIndFromNumFeat]);

%initialize variables where results will be stored
motionStatsSample = repmat(struct('probMotionType',[],'motionChar',[]),...
    maxNumRep,numCells);
msStatsSample = repmat(struct('statsGeneral',[],'motionTypeFrac',[],...
    'probMSIfTypeFull',[],'probMSIfType12',[],'numEventsPerCatFull',[],...
    'motionRedTypeFrac',[],'probMSIfRedTypeFull',[],'probMSIfRedType12',[]),...
    maxNumRep,numCells);

%go over these cells
for iCell = 1 : numCells
    
    %put the original cell statistics as the first entry
    motionStatsSample(1,iCell) = motionStatsInd(iCell);
    msStatsSample(1,iCell) = msStatsInd(iCell);
    
    %remove tracks shorter than minTrackLen
    tracksTmp = trackData(iCell).tracks;
    criteria.lifeTime.min = minTrackLen;
    indxKeep = chooseTracks(tracksTmp,criteria);
    tracksTmp = tracksTmp(indxKeep);
    diffAnalysisTmp = diffAnalysisInd(iCell).results(indxKeep);
    
    %get number of remaining tracks in this cell
    numTracks = length(tracksTmp);
    
    %generate (number of repetitions - 1) bootstrap samples for this cell
    %to fill the rest of its entries
    for iBoot = 2 : max([numRepIndFromNumFeat(iCell) numRepIndFromNumFeatSquared(iCell)])
        
        %select a random subset of tracks with replacement
        randIndx = randsample((1:numTracks),numTracks,'true');
        tracksSample = tracksTmp(randIndx);
        diffAnalysisSample = diffAnalysisTmp(randIndx);
        
        %get motion statistics
        [motionStatsSample(iBoot,iCell).probMotionType,motionStatsSample(iBoot,iCell).motionChar] = ...
            summarizeDiffAnRes(tracksSample,minTrackLen,probDim,...
            diffAnalysisSample,extractType);
        
        %get merge/split statistics
        [msStatsSample(iBoot,iCell).statsGeneral,msStatsSample(iBoot,iCell).motionTypeFrac,...
            msStatsSample(iBoot,iCell).probMSIfTypeFull,msStatsSample(iBoot,iCell).probMSIfType12,...
            msStatsSample(iBoot,iCell).numEventsPerCatFull,msStatsSample(iBoot,iCell).motionRedTypeFrac,...
            msStatsSample(iBoot,iCell).probMSIfRedTypeFull,msStatsSample(iBoot,iCell).probMSIfRedType12] = ...
            calcStatsMS(tracksSample,minTrackLen,probDim,diffAnalysisSample,removePotArtifacts);
        
    end
    
end

%% Generate weighted distribution of merging and splitting statistics using the bootstrap samples

%initialize variable
probMSIfType1Distr = []; %m/s probability calculated using Method 1
probMSIfType2Distr = []; %m/s probability calculated using Method 2

%go over all cells
for iCell = 1 : numCells
    
    %get this cell's entries
    msStatsTmp = msStatsSample(1:numRepIndFromNumFeatSquared(iCell),iCell); %#ok<NASGU>
    probMSIfType1Ind = (catStruct(2,'msStatsTmp.probMSIfType12(:,1,3)'))';
    probMSIfType2Ind = (catStruct(2,'msStatsTmp.probMSIfType12(:,2,3)'))';
    
    %add its original contribution to the weighted distribution
    probMSIfType1Distr = [probMSIfType1Distr; probMSIfType1Ind]; %#ok<AGROW>
    probMSIfType2Distr = [probMSIfType2Distr; probMSIfType2Ind]; %#ok<AGROW>
    
end

%% Generate weighted distribution of motion type statistics using the bootstrap samples

%initialize variable
probMotionTypeDistr = []; %motion type probabilities
linearDiffCoefDistr = []; %diffusion coefficient of linear particles
linearConfWidthDistr = []; %width of linear confinement areas
brownianDiffCoefDistr = []; %diffusion coefficient of brownian particles
confinedDiffCoefDistr = []; %diffusion cofficient of confined particles
confinedConfWidthDistr = []; %size of confinement area for confined particles

%go over all cells
for iCell = 1 : numCells
    
    %get this cell's entries
    motionStatsTmp = motionStatsSample(1:numRepIndFromNumFeat(iCell),iCell); %#ok<NASGU>
    probMotionTypeInd = (catStruct(2,'motionStatsTmp.probMotionType(:,1)'))';
    linearDiffCoefInd = (catStruct(2,'motionStatsTmp.motionChar.linear.all.meanStd(:,1)'))';
    linearConfWidthInd = (catStruct(2,'motionStatsTmp.motionChar.linear.all.meanStd(:,2)'))';
    brownianDiffCoefInd = (catStruct(2,'motionStatsTmp.motionChar.notLinear.brownian.meanStd(:,1)'))';
    confinedDiffCoefInd = (catStruct(2,'motionStatsTmp.motionChar.notLinear.confined.meanStd(:,1)'))';
    confinedConfWidthInd = [];
    confinedConfWidthInd(:,:,1) = (catStruct(2,'motionStatsTmp.motionChar.notLinear.confined.meanStd(:,2)'))';
    confinedConfWidthInd(:,:,2) = (catStruct(2,'motionStatsTmp.motionChar.notLinear.confined.meanStd(:,3)'))';
    confinedConfWidthInd = mean(confinedConfWidthInd,3);
    
    probMotionTypeDistr = [probMotionTypeDistr; probMotionTypeInd]; %#ok<AGROW>
    linearDiffCoefDistr = [linearDiffCoefDistr; linearDiffCoefInd]; %#ok<AGROW>
    linearConfWidthDistr = [linearConfWidthDistr; linearConfWidthInd]; %#ok<AGROW>
    brownianDiffCoefDistr = [brownianDiffCoefDistr; brownianDiffCoefInd]; %#ok<AGROW>
    confinedDiffCoefDistr = [confinedDiffCoefDistr; confinedDiffCoefInd]; %#ok<AGROW>
    confinedConfWidthDistr = [confinedConfWidthDistr; confinedConfWidthInd]; %#ok<AGROW>
    
end

%% Output variables

%full analysis results
analysisResFull.all.original = struct('diffusion',diffAnalysisAll,...
    'msStats',msStatsAll,'motionStats',motionStatsAll);
analysisResFull.all.bootstrap = struct('msStats',msStatsAllSample,...
    'motionStats',motionStatsAllSample);
analysisResFull.ind.original = struct('diffusion',diffAnalysisInd,...
    'msStats',msStatsInd,'motionStats',motionStatsInd);
analysisResFull.ind.bootstrap = struct('msStats',msStatsSample,...
    'motionStats',motionStatsSample,...
    'numRepFromNumFeat',numRepIndFromNumFeat,...
    'numRepFromNumFeatSquared',numRepIndFromNumFeatSquared);

%results of particular interest

%motion type probabilities
analysisResPart.probMotionType.all.original  = motionStatsAll.probMotionType(:,1)';
analysisResPart.probMotionType.all.bootstrap = (catStruct(2,'motionStatsAllSample.probMotionType(:,1)'))';
analysisResPart.probMotionType.ind.original  = (catStruct(2,'motionStatsInd.probMotionType(:,1)'))';
analysisResPart.probMotionType.ind.bootstrap = probMotionTypeDistr;

%probability of linear motion
analysisResPart.probLin.all.original  = sum(analysisResPart.probMotionType.all.original([1 2 3 4 7]));
analysisResPart.probLin.all.bootstrap = sum(analysisResPart.probMotionType.all.bootstrap(:,[1 2 3 4 7]),2);
analysisResPart.probLin.ind.original  = sum(analysisResPart.probMotionType.ind.original(:,[1 2 3 4 7]),2);
analysisResPart.probLin.ind.bootstrap = sum(analysisResPart.probMotionType.ind.bootstrap(:,[1 2 3 4 7]),2);

%diffusion coefficient of linear particles
analysisResPart.linearDiffCoef.all.original  = motionStatsAll.motionChar.linear.all.meanStd(:,1)';
analysisResPart.linearDiffCoef.all.bootstrap = (catStruct(2,'motionStatsAllSample.motionChar.linear.all.meanStd(:,1)'))';
analysisResPart.linearDiffCoef.ind.original  = (catStruct(2,'motionStatsInd.motionChar.linear.all.meanStd(:,1)'))';
analysisResPart.linearDiffCoef.ind.bootstrap = linearDiffCoefDistr;

%confinement width of linear particles
analysisResPart.linearConfHalfWidth.all.original  = motionStatsAll.motionChar.linear.all.meanStd(:,2)';
analysisResPart.linearConfHalfWidth.all.bootstrap = (catStruct(2,'motionStatsAllSample.motionChar.linear.all.meanStd(:,2)'))';
analysisResPart.linearConfHalfWidth.ind.original  = (catStruct(2,'motionStatsInd.motionChar.linear.all.meanStd(:,2)'))';
analysisResPart.linearConfHalfWidth.ind.bootstrap = linearConfWidthDistr;

%diffusion coefficient of Brownian particles
analysisResPart.brownianDiffCoef.all.original  = motionStatsAll.motionChar.notLinear.brownian.meanStd(:,1)';
analysisResPart.brownianDiffCoef.all.bootstrap = (catStruct(2,'motionStatsAllSample.motionChar.notLinear.brownian.meanStd(:,1)'))';
analysisResPart.brownianDiffCoef.ind.original  = (catStruct(2,'motionStatsInd.motionChar.notLinear.brownian.meanStd(:,1)'))';
analysisResPart.brownianDiffCoef.ind.bootstrap = brownianDiffCoefDistr;

%diffusion coefficient of confined particles
analysisResPart.confinedDiffCoef.all.original  = motionStatsAll.motionChar.notLinear.confined.meanStd(:,1)';
analysisResPart.confinedDiffCoef.all.bootstrap = (catStruct(2,'motionStatsAllSample.motionChar.notLinear.confined.meanStd(:,1)'))';
analysisResPart.confinedDiffCoef.ind.original  = (catStruct(2,'motionStatsInd.motionChar.notLinear.confined.meanStd(:,1)'))';
analysisResPart.confinedDiffCoef.ind.bootstrap = confinedDiffCoefDistr;

%confinement width of confined particles
confinedConfWidthInd = [];
confinedConfWidthInd(:,:,1) = motionStatsAll.motionChar.notLinear.confined.meanStd(:,2)';
confinedConfWidthInd(:,:,2) = motionStatsAll.motionChar.notLinear.confined.meanStd(:,3)';
confinedConfWidthInd = mean(confinedConfWidthInd,3);
analysisResPart.confinedConfHalfWidth.all.original  = confinedConfWidthInd;
confinedConfWidthInd = [];
confinedConfWidthInd(:,:,1) = (catStruct(2,'motionStatsAllSample.motionChar.notLinear.confined.meanStd(:,2)'))';
confinedConfWidthInd(:,:,2) = (catStruct(2,'motionStatsAllSample.motionChar.notLinear.confined.meanStd(:,3)'))';
confinedConfWidthInd = mean(confinedConfWidthInd,3);
analysisResPart.confinedConfHalfWidth.all.bootstrap = confinedConfWidthInd;
confinedConfWidthInd = [];
confinedConfWidthInd(:,:,1) = (catStruct(2,'motionStatsInd.motionChar.notLinear.confined.meanStd(:,2)'))';
confinedConfWidthInd(:,:,2) = (catStruct(2,'motionStatsInd.motionChar.notLinear.confined.meanStd(:,3)'))';
confinedConfWidthInd = mean(confinedConfWidthInd,3);
analysisResPart.confinedConfHalfWidth.ind.original  = confinedConfWidthInd;
analysisResPart.confinedConfHalfWidth.ind.bootstrap = confinedConfWidthDistr;

%merge/split probabilities - 1
analysisResPart.probMSIfType1.all.original = msStatsAll.probMSIfType12(:,1,3)';
analysisResPart.probMSIfType1.all.bootstrap = (catStruct(2,'msStatsAllSample.probMSIfType12(:,1,3)'))';
analysisResPart.probMSIfType1.ind.original = (catStruct(2,'msStatsInd.probMSIfType12(:,1,3)'))';
analysisResPart.probMSIfType1.ind.bootstrap = probMSIfType1Distr;

%merge/split probabilities - 2
analysisResPart.probMSIfType2.all.original = msStatsAll.probMSIfType12(:,2,3)';
analysisResPart.probMSIfType2.all.bootstrap = (catStruct(2,'msStatsAllSample.probMSIfType12(:,2,3)'))';
analysisResPart.probMSIfType2.ind.original = (catStruct(2,'msStatsInd.probMSIfType12(:,2,3)'))';
analysisResPart.probMSIfType2.ind.bootstrap = probMSIfType2Distr;

