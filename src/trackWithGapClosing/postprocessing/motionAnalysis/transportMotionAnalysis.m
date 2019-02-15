function [motionAnalysisWithErrBar,motionAnalysisNoErrBar,errFlag] = ...
    transportMotionAnalysis(tracks,refCoord,pixelSize,timeBetweenFrames,...
    extractType,probDim,diffAnalysisAlphaValues,verbose,minLength)
%TRANSPORTMOTIONANALYSIS extracts the properties of motor-driven motion
%
%SYNOPSIS [motionAnalysisWithErrBar,motionAnalysisNoErrBar,errFlag] = ...
%    transportMotionAnalysis(tracks,refCoord,pixelSize,timeBetweenFrames,...
%    extractType,probDim,diffAnalysisAlphaValues,verbose)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits.
%                     Cose assumes that coordinates are in pixels.
%       refCoord    : Coordinates of reference point to assign
%                     positive and negative directions, in pixels.
%       pixelSize   : Pixel size, in whatever units user wants.
%                     Optional. Default: 1.
%       timeBetweenFrames: Time between frames, in whatever units user
%                     wants.
%                     Optional. Default: 1.
%       extractType : 1 - Analyze every track segment separately.
%                     2 - Extract from each compound track the longest
%                         trajectory to use in analysis - NOT IMPLEMENTED
%                         YET.
%                     Variable irrelevant if tracks are input as a matrix.
%                     Optional. Default: 1.
%       probDim     : Problem dimensionality. Optional. Default: 2.
%       diffAnalysisAlphaValues: Row vector with 2 entries. First entry is
%                     the alpha-value for MSS analysis (can take the values
%                     0.2, 0.1, 0.05 and 0.01; see help of trackMSSAnalysis
%                     for most up-to-date values allowed). Second entry is
%                     the alpha-value for asymmetry determination (can take
%                     the values 0.2, 0.1, 0.05 and 0.01; see help of
%                     asymDeterm2D3D for most up-to-date values allowed).
%                     Optional. Default: [0.05 0.2]. If only one value is
%                     entered, it is taken as the alpha-value for MSS
%                     analysis, while the alpha-value for asymmetry
%                     analysis is given the default value.
%       verbose     : 0 to make no plots, 1 to show plots of final results,
%                     2 to show also plots of intermediate results.
%                     Optional. Default: 1.
%       minLength   : Minimum length of tracks to analyze. 
%                     Optional. Default: 1.
%
%OUTPUT motionAnalysisWithErrBar: Structure storing the results of
%                     analyzing the motion taking positional uncertainties
%                     into account. Contains 3 fields:
%             .awayFromReferencePoint: Contains the sub-fields
%                  .instantSpeed: Instantaneous speed.
%                  .runTime     : The time of each run.
%                  .runLength   : The length of each run.
%                  .averageSpeed: runLength / runTime.
%                     Note that "runTime", "runLength" and "averageSpeed"
%                     have the same number of values, while "instantSpeed"
%                     is different.
%             .towardReferencePoint  : Contains the same sub-fields as above.
%             .pause                 : Contains 1 sub-field
%                  .pauseTime"  : Duration of each pause.
%                     Each of the above sub-fields is a structure array
%                     with 10 entries, storing the results for
%                     groups of trajectories. The groups are based on
%                     motion type, and are the following:
%                       (1) Linear + super-diffusion
%                       (2) Linear + normal-diffusion
%                       (3) Linear + sub-diffusion
%                       (4) Linear + undetermined diffusion
%                       (5) Not linear + super-diffusion
%                       (6) Not linear + normal-diffusion
%                       (7) Not linear + sub-diffusion
%                       (8) All linear, i.e. 1+2+3+4
%                       (9) All not linear, i.e. 5+6+7
%                       (10) All trajectories, i.e. 1+2+3+4+5+6+7.
%                     Each entry contains 2 sub-sub-fields:
%                       .values: Listing all the values.
%                       .stats : Showing the mean, median, std, min and max
%                                of the distribution of values.
%       motionAnalysisNoErrBar: Structure storing the results of analyzing
%                     the motion ignoring positional uncertainties.
%                     Contains the fields "awayFromReferencePoint" and
%                     "towardReferencePoint". Inside them, there is only
%                     the sub-fields "instantSpeed" and "runTime", which
%                     are identical in organization to the fields with the
%                     same names in "motionAnalysisWithErrBar".
%
%REMARK Output speeds and times will be in whatever units user enters for
%pixel size and time between frames. If no units are entered, then the
%output will be in pixels and frames.
%
%Khuloud Jaqaman, December 2010

%% Output

motionAnalysisWithErrBar = [];
motionAnalysisNoErrBar = [];
errFlag = 0;

%% Input

%check whether tracks were input
if nargin < 2
    disp('--transportMotionAnalysis: Mandatory first two input arguments missing!');
    errFlag = 1;
    return
end

if nargin < 3 || isempty(pixelSize)
    pixelSize = 1;
end

if nargin < 4 || isempty(timeBetweenFrames)
    timeBetweenFrames = 1;
end

if nargin < 5 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--transportMotionAnalysis: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

if nargin < 6 || isempty(probDim)
    probDim = 2;
end

if nargin < 7 || isempty(diffAnalysisAlphaValues)
    diffAnalysisAlphaValues = [0.05 0.2];
elseif length(diffAnalysisAlphaValues) == 1
    diffAnalysisAlphaValues = [diffAnalysisAlphaValues 0.2];
end

if nargin < 8 || isempty(verbose)
    verbose = 1;
end

if nargin < 9 || isempty(minLength)
    minLength = 1;
end

if errFlag
    disp('--transportMotionAnalysis: Please fix input variables');
    return
end

%% Track extraction for analysis

%keep for analysis only tracks of minimum length
criteria.lifeTime.min = minLength;
indxKeep = chooseTracks(tracks,criteria);
tracks = tracks(indxKeep);

%store input tracks in a new variable
tracksInput = tracks;

%extract segments for analysis if tracks were input as a structure that
%might contain merges and splits
%the point is to reduce compound tracks that contain merges and splits into
%simple separate tracks
%thus this step is not necessary if the tracks were input as a matrix,
%which by definition does not contain unresolved compound tracks.
if isstruct(tracks)
    
    %get number of input tracks from structure
    numInputTracks = length(tracksInput);
    
    clear tracks
    
    switch extractType
        
        case 1 %retrieve every track segment separately
            
            [tracks,dummy,compTrackStartRow,numSegments] = ...
                convStruct2MatIgnoreMS(tracksInput);
            
        case 2 %make the longest track possible, given all the merges and splits
            
            disp('Sorry - not implement yet!')
            errFlag = 1;
            return
            
    end
    
else
    
    %get number of input tracks from matrix
    numInputTracks = size(tracksInput,1);
    
    %indicate rows where tracks start (trivial in this case)
    compTrackStartRow = (1 : numInputTracks)';
    
    %indicate number of segments in each track (1 for all tracks)
    numSegments = ones(numInputTracks,1);
    
end

%get number of track segments to be analyzed
numTrackSegments = size(tracks,1);

%% Diffusion analysis to remove tracks that are most probably not motor-driven

%perform diffusion analysis
[diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,extractType,...
    probDim,1,diffAnalysisAlphaValues,0,0);

%get track classification
trackClass = vertcat(diffAnalysisRes.classification);

%find indices of tracks to keep - those are tracks that are either
%asymmetric or exhibit directed motion
indxKeep = find(trackClass(:,1)==1 | trackClass(:,2)==3);
indxRemove = setdiff(1:numTrackSegments,indxKeep)';

%separate tracks to keep from those to remove
tracksRemove = tracks(indxRemove,:);
tracks = tracks(indxKeep,:);

if isempty(tracks)
    disp('No tracks to analyze. Exiting ...')
    return
end

%retain only the diffusion analysis results of the tracks to keep
diffAnalysisResRemove = diffAnalysisRes(indxRemove);
diffAnalysisRes = diffAnalysisRes(indxKeep);
% trackClass = trackClass(indxKeep,:);

%show tracks kept and tracks removed for visual assessment
if verbose >= 2
    if ~isempty(tracksRemove)
        %         figure('Name','Removed tracks')
        plotTracksDiffAnalysis2D(tracksRemove,diffAnalysisResRemove,[],1,[],0,0);
        %         plotTracks2D(tracksRemove,[],'2',[],0,1,[],0,0);
    else
        disp('No tracks removed');
    end
    %     figure('Name','Retained tracks to be analyzed')
    %     plotTracksDiffAnalysis2D(tracks,diffAnalysisRes,[],1,[],0,0);
    plotTracks2D(tracks,[],'1',[],1,1,[],0,0);
    hold on
    plot(refCoord(1)*[1 1],refCoord(2)*[1 1],'kx','LineWidth',2,'MarkerSize',10)
end

%% Simple analysis

%perform my simple run time and velocity analysis
%importantly, get from it the positions along the direction of motion
[runTimePos,runTimeNeg,dispPos,dispNeg,runTimePerp,dispPerp,posAlongDir] = ...
    runLengthAnalysis(tracks,refCoord,5,diffAnalysisRes);

%find number of tracks in each motion category and determine types that
%have at least one track
numTracksPerType = zeros(7,1);
for iType = 1 : 7
    numTracksPerType(iType) = length(posAlongDir(iType).values);
end
indxGoodType = find(numTracksPerType>0);

%% Analysis that takes error bars into account

%convert position-along-direction-of-motion time series to run Jonas's
%trajectory analysis
posAlongDirJ = repmat(struct('values',[]),7,1);
for iType = indxGoodType'
    posAlongDirJ(iType).values = convertTrajectoryData(posAlongDir(iType).values);
end

%run Jonas's trajectory analysis
statsJ = repmat(struct('values',[]),7,1);
for iType = indxGoodType'
    statsJ(iType).values = trajectoryAnalysis(posAlongDirJ(iType).values,struct('verbose',0,'saveTxt',0));
end

%run the comparison code in order to extract speed and time distributions
[globCompare,globStats]=testTrajectories([statsJ(1).values statsJ(2).values ...
    statsJ(3).values statsJ(4).values statsJ(5).values statsJ(6).values ...
    statsJ(7).values],0);
statsJExpanded(indxGoodType) = globStats;

%get the run lengths in both directions
for iType = indxGoodType'
    
    runLengthPosTmp = [];
    runLengthNegTmp = [];
    numReversalsTmp = [];
    for iTrack = 1 : numTracksPerType(iType)
        
        %get dataListGroup from output of trajectoryAnalysis
        dataListG = statsJ(iType).values.individualStatistics(iTrack).dataListGroup;
        
        %get intervals away from reference point
        [dummy,posGroups] = trajectoryAnalysisMainCalcStatsFindGroups(dataListG,1);
        
        %calculate run lengths in these intervals
        if ~isempty(dummy)
            tmptmp = [];
            for iInt = 1 : size(posGroups,1)
                tmptmp = [tmptmp; sum(dataListG(posGroups(iInt,1):posGroups(iInt,2),9))];
            end
            runLengthPosTmp = [runLengthPosTmp; tmptmp];
        end
        
        %get intervals toward reference point
        [dummy,negGroups] = trajectoryAnalysisMainCalcStatsFindGroups(dataListG,2);
        
        %calculate run lengths in these intervals
        if ~isempty(dummy)
            tmptmp = [];
            for iInt = 1 : size(negGroups,1)
                tmptmp = [tmptmp; sum(dataListG(negGroups(iInt,1):negGroups(iInt,2),9))];
            end
            runLengthNegTmp = [runLengthNegTmp; -tmptmp];
        end
        
        %calculate number of transitions in this trajectory
        trajHist = dataListG(:,3);
        trajHist = trajHist(trajHist==1|trajHist==2);
        trajHistDiff = diff(trajHist);
        numReversalsTmp = [numReversalsTmp; length(find(trajHistDiff~=0))];
        
    end
    
    %save data from all tracks belonging to this motion type
    runLengthPos(iType).values = runLengthPosTmp;
    runLengthNeg(iType).values = runLengthNegTmp;
    numReversals(iType).values = numReversalsTmp;
    
    
    %calculate the average speed for each run
    averageSpeedPos(iType).values = runLengthPosTmp./statsJExpanded(iType).growthTimes;
    averageSpeedNeg(iType).values = runLengthNegTmp./statsJExpanded(iType).shrinkageTimes;
    
end

%% Collect results for output

%analysis results that ignore error bars
awayFromRefPoint = struct('instantSpeed',dispPos,'runTime',runTimePos);
towardRefPoint = struct('instantSpeed',dispNeg,'runTime',runTimeNeg);
motionAnalysisNoErrBar = struct('awayFromRefPoint',awayFromRefPoint,...
    'towardRefPoint',towardRefPoint);

%analysis results that use error bars
awayFromRefPoint = struct('instantSpeed',runTimePos,'runTime',runTimePos,...
    'runLength',runTimePos,'averageSpeed',runTimePos);
towardRefPoint = awayFromRefPoint;
pauseInfo = struct('pauseTime',runTimePos);
% reversalInfo = struct('numPerTraj',[]);
for iType = indxGoodType'
    awayFromRefPoint.instantSpeed(iType).values = statsJExpanded(iType).growthSpeeds/60;
    awayFromRefPoint.runTime(iType).values = statsJExpanded(iType).growthTimes;
    towardRefPoint.instantSpeed(iType).values = statsJExpanded(iType).shrinkageSpeeds/60;
    towardRefPoint.runTime(iType).values = statsJExpanded(iType).shrinkageTimes;
    tmp = statsJ(iType).values.overallDistribution.pauseTime;
    if ~isempty(tmp)
        pauseInfo.pauseTime(iType).values = tmp(:,1);
    end
end
awayFromRefPoint.runLength(indxGoodType) = runLengthPos(indxGoodType);
awayFromRefPoint.averageSpeed(indxGoodType) = averageSpeedPos(indxGoodType);
towardRefPoint.runLength(indxGoodType) = runLengthNeg(indxGoodType);
towardRefPoint.averageSpeed(indxGoodType) = averageSpeedNeg(indxGoodType);

motionAnalysisWithErrBar.awayFromRefPoint = awayFromRefPoint;
motionAnalysisWithErrBar.towardRefPoint = towardRefPoint;
motionAnalysisWithErrBar.pause = pauseInfo;

%each motion characteristic so far is divided into 7 motion categories:
%   (1) linear + super-diffusive
%   (2) linear + normal-diffusive
%   (3) linear + sub-diffusive
%   (4) linear + undetermined
%   (5) nonlinear + super-diffusive
%   (6) nonlinear + normal-diffusive
%   (7) nonlinear + sub-diffusive
%now merge categories to get the following:
%   (8) all linear
%   (9) all nonlinear
%   (10) all linear + all nonlinear = everything
varAwayToward = {'awayFromRefPoint','towardRefPoint'};
varMotionChar = {'instantSpeed','runTime','runLength','averageSpeed'};
%analysis results that ignore error bars
for iVarAT = 1 : 2
    for iVarMC = 1 : 2
        eval(['tmp1 = vertcat(motionAnalysisNoErrBar.' ...
            varAwayToward{iVarAT} '.' varMotionChar{iVarMC} '(1:4).values);'])
        eval(['tmp2 = vertcat(motionAnalysisNoErrBar.' ...
            varAwayToward{iVarAT} '.' varMotionChar{iVarMC} '(5:7).values);'])
        eval(['motionAnalysisNoErrBar.' varAwayToward{iVarAT} '.' ...
            varMotionChar{iVarMC} '(8).values = tmp1;'])
        eval(['motionAnalysisNoErrBar.' varAwayToward{iVarAT} '.' ...
            varMotionChar{iVarMC} '(9).values = tmp2;'])
        eval(['motionAnalysisNoErrBar.' varAwayToward{iVarAT} '.' ...
            varMotionChar{iVarMC} '(10).values = [tmp1; tmp2];'])
    end
end
%analysis results that use error bars
for iVarAT = 1 : 2
    for iVarMC = 1 : 4
        eval(['tmp1 = vertcat(motionAnalysisWithErrBar.' ...
            varAwayToward{iVarAT} '.' varMotionChar{iVarMC} '(1:4).values);'])
        eval(['tmp2 = vertcat(motionAnalysisWithErrBar.' ...
            varAwayToward{iVarAT} '.' varMotionChar{iVarMC} '(5:7).values);'])
        eval(['motionAnalysisWithErrBar.' varAwayToward{iVarAT} '.' ...
            varMotionChar{iVarMC} '(8).values = tmp1;'])
        eval(['motionAnalysisWithErrBar.' varAwayToward{iVarAT} '.' ...
            varMotionChar{iVarMC} '(9).values = tmp2;'])
        eval(['motionAnalysisWithErrBar.' varAwayToward{iVarAT} '.' ...
            varMotionChar{iVarMC} '(10).values = [tmp1; tmp2];'])
    end
end
tmp1 = vertcat(motionAnalysisWithErrBar.pause.pauseTime(1:4).values);
tmp2 = vertcat(motionAnalysisWithErrBar.pause.pauseTime(5:7).values);
motionAnalysisWithErrBar.pause.pauseTime(8).values = tmp1;
motionAnalysisWithErrBar.pause.pauseTime(9).values = tmp2;
motionAnalysisWithErrBar.pause.pauseTime(10).values = [tmp1;tmp2];

%calculate for each distribution the mean, median, std, min and max
%also convert to user-specified units
unitConversion = [pixelSize/timeBetweenFrames timeBetweenFrames ...
    pixelSize pixelSize/timeBetweenFrames]';
%analysis results that ignore error bars
for iVarAT = 1 : 2
    for iVarMC = 1 : 2
        for iType = 1 : 10
            eval(['tmp = motionAnalysisNoErrBar.' varAwayToward{iVarAT} ...
                '.' varMotionChar{iVarMC} '(iType).values;'])
            if ~isempty(tmp)
                tmp = tmp*unitConversion(iVarMC);
                tmpStats = [mean(tmp) median(tmp) std(tmp) min(tmp) max(tmp)];
                eval(['motionAnalysisNoErrBar.' varAwayToward{iVarAT} ...
                    '.' varMotionChar{iVarMC} '(iType).values = tmp;'])
                eval(['motionAnalysisNoErrBar.' varAwayToward{iVarAT} ...
                    '.' varMotionChar{iVarMC} '(iType).stats = tmpStats;'])
            else
                eval(['motionAnalysisNoErrBar.' varAwayToward{iVarAT} ...
                    '.' varMotionChar{iVarMC} '(iType).stats = [];'])
            end
        end
    end
end
%analysis results that use error bars
for iVarAT = 1 : 2
    for iVarMC = 1 : 4
        for iType = 1 : 10
            eval(['tmp = motionAnalysisWithErrBar.' varAwayToward{iVarAT} ...
                '.' varMotionChar{iVarMC} '(iType).values;'])
            if ~isempty(tmp)
                tmp = tmp*unitConversion(iVarMC);
                tmpStats = [mean(tmp) median(tmp) std(tmp) min(tmp) max(tmp)];
                eval(['motionAnalysisWithErrBar.' varAwayToward{iVarAT} ...
                    '.' varMotionChar{iVarMC} '(iType).values = tmp;'])
                eval(['motionAnalysisWithErrBar.' varAwayToward{iVarAT} ...
                    '.' varMotionChar{iVarMC} '(iType).stats = tmpStats;'])
            else
                eval(['motionAnalysisWithErrBar.' varAwayToward{iVarAT} ...
                    '.' varMotionChar{iVarMC} '(iType).stats = [];'])
            end
        end
    end
end
for iType = 1 : 10
    tmp = motionAnalysisWithErrBar.pause.pauseTime(iType).values;
    if ~isempty(tmp)
        tmp = tmp*timeBetweenFrames;
        tmpStats = [mean(tmp) median(tmp) std(tmp) min(tmp) max(tmp)];
        motionAnalysisWithErrBar.pause.pauseTime(iType).values = tmp;
        motionAnalysisWithErrBar.pause.pauseTime(iType).stats = tmpStats;
    else
        motionAnalysisWithErrBar.pause.pauseTime(iType).stats = [];
    end
end

%% Plot results

if verbose >= 1
    plotTransportMotionAnalysisRes(motionAnalysisWithErrBar,10);
end

%% ~~~ the end ~~~

