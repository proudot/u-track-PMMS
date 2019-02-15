function [msdAnalysis,mssAnalysis,projAnalysis,errFlag] = trackDiffAnalysisMod(tracks,...
    analysisType,extractType,probDim)
%TRACKDIFFANALYSISMOD performs various diffusion analyses on input tracks
%
%SYNOPSIS [msdAnalysis,mssAnalysis,projAnalysis,errFlag] = trackDiffAnalysisMod(tracks,...
%    analysisType,extractType,probDim)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits.
%       analysisType: A vector containing any one, or combination, of the
%                     following options:
%                     1 - Mean square displacement analysis of tracks in a
%                         manner similar to the Huet paper (BJ 2006),
%                         although without windowing for now.
%                     2 - Moments scaling spectrum analysis of tracks using
%                         the method in Ewers (PNAS 2005) and Ferrari
%                         (Physica D 2001).
%                     3 - Moments scaling spectrum analysis on the
%                         PROJECTION of linear tracks onto their axis of
%                         preferred direction.
%                     Optional. Default: [1 2].
%       extractType : 1 - Analyze every track segment separately.
%                     2 - Extract from each compound track the longest
%                         trajectory to use in analysis - NOT IMPLEMENTED
%                         YET.
%                     Variable irrelevant if tracks are input as a matrix.
%                     Optional. Default: 1.
%       probDim     : Problem dimensionality. Optional. Default: 2.
%
%OUTPUT msdAnalysis : Results of Analysis 1.
%       mssAnalysis : Results of Analysis 2.
%       projAnalysis: Results of Analysis 3.
%       errFlag     : 0 if executed normally, 1 otherwise.
%
%REMARKS AnalysisType 1 cannot be applied to 1D data.
%
%Khuloud Jaqaman, February 2008

%% Output

msdAnalysis = [];
mssAnalysis = [];
projAnalysis = [];
errFlag = 0;

%% Input

%check whether tracks were input
if nargin < 1
    disp('--trackDiffAnalysisMod: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(analysisType)
    analysisType = [1 2];
else
    if ~all(analysisType == 1 | analysisType == 2 | analysisType == 3)
        disp('--trackDiffAnalysisMod: Variable analysisType should include only 1, 2 and 3.');
        errFlag = 1;
    end
    if any(analysisType == 3) && ~any(analysisType == 1 | analysisType == 2)
        analysisType(end+1) = 1;
    end
end

if nargin < 3 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--trackDiffAnalysisMod: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end
    
if nargin < 4 || isempty(probDim)
    probDim = 2;
end

if errFlag
    disp('--trackDiffAnalysisMod: Please fix input variables');
    return
end

%% track extraction for analysis

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

%% analysis

%Mean square displacement analysis on tracks (Huet et al, BJ 2006)
if any(analysisType == 1)

    %do the analysis
    [trackClass,diffCoef,confParam,asymParam] = trackMSDAnalysis(...
        tracks,probDim);

    %save the analysis
    msdAnalysis = repmat(struct('classification',[],'diffCoef',[],...
        'confParam',[],'asymParam',[]),numInputTracks,1);
    for iTrack = 1 : numInputTracks
        msdAnalysis(iTrack).classification = trackClass(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        msdAnalysis(iTrack).diffCoef = diffCoef(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        msdAnalysis(iTrack).confParam = confParam(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        msdAnalysis(iTrack).asymParam = asymParam(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
    end

end

%Moment scaling spectrum analysis on tracks (Ewers et al, PNAS 2005 &
%Ferrari et al, Physica D 2001)
if any(analysisType == 2)

    %do the analysis
    momentOrders = 0 : 6;
    [trackClass,mssSlope,genDiffCoef,scalingPower] = trackMSSAnalysis(...
        tracks,probDim,momentOrders);

    %save the analysis
    mssAnalysis = repmat(struct('classification',[],'momentScalingPowerSlope',[],...
        'momentOrders',momentOrders,'genDiffCoef',[],'momentScalingPower',[]),...
        numInputTracks,1);
    for iTrack = 1 : numInputTracks
        mssAnalysis(iTrack).classification = trackClass(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        mssAnalysis(iTrack).momentScalingPowerSlope = mssSlope(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        mssAnalysis(iTrack).genDiffCoef = genDiffCoef(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1,:);
        mssAnalysis(iTrack).momentScalingPower = scalingPower(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1,:);
    end

end

%Moment scaling spectrum analysis on projection of linear tracks
if any(analysisType == 3)
    
    %get linear track segments
    if any(analysisType == 1)
        segmentClass = vertcat(msdAnalysis.classification);
    elseif any(analysisType == 2)
        segmentClass = vertcat(mssAnalysis.classification);
    else
        disp('--trackDiffAnalysisMod: Please allow classification of tracks in 2D/3D before classifying the linear ones in 1D.');
        errFlag = 1;
        return
    end
    indxLin = find(segmentClass == 3);

    %project positions of linear tracks onto direction of motion
    numCol = size(tracks,2) / 8;
    tracksLin = NaN(length(indxLin),numCol*8);
    iLin = 0;
    for iTrack = indxLin'

        %get the positions in this track and their standard deviations
        trackCoordX = tracks(iTrack,1:8:end)';
        deltaCoordX = tracks(iTrack,5:8:end)';
        trackCoordY = tracks(iTrack,2:8:end)';
        deltaCoordY = tracks(iTrack,6:8:end)';
        trackCoordZ = tracks(iTrack,3:8:end)';
        deltaCoordZ = tracks(iTrack,7:8:end)';
        trackCoord = [trackCoordX trackCoordY trackCoordZ];
        deltaCoord = [deltaCoordX deltaCoordY deltaCoordZ];
        trackCoord = trackCoord(:,1:probDim);
        deltaCoord = deltaCoord(:,1:probDim);

        %project onto direction of motion
        [posAlongDir,deltaPosAlongDir] = projectCoordOntoDir(trackCoord,...
            deltaCoord,[],[]);

        %construct matrix of linear tracks with projected positions
        trackCoord2 = [posAlongDir zeros(numCol,3) deltaPosAlongDir zeros(numCol,3)]';
        trackCoord2 = trackCoord2(:)';
        iLin = iLin + 1;
        tracksLin(iLin,:) = trackCoord2;

    end

    %do the analysis
    momentOrders = 0 : 6;
    [trackClassT,mssSlopeT,genDiffCoefT,scalingPowerT] = trackMSSAnalysis(...
        tracksLin,1,momentOrders);
    
    %since not all track segments are linear, put analysis results in their
    %proper place among all track segment
    trackClass = NaN(numTrackSegments,1);
    mssSlope = trackClass;
    genDiffCoef = NaN(numTrackSegments,length(momentOrders));
    scalingPower = genDiffCoef;
    trackClass(indxLin) = trackClassT;
    mssSlope(indxLin) = mssSlopeT;
    genDiffCoef(indxLin,:) = genDiffCoefT;
    scalingPower(indxLin,:) = scalingPowerT;

    %save the analysis
    projAnalysis = repmat(struct('classification',[],'momentScalingPowerSlope',[],...
        'momentOrders',momentOrders,'genDiffCoef',[],'momentScalingPower',[]),...
        numInputTracks,1);
    for iTrack = 1 : numInputTracks
        projAnalysis(iTrack).classification = trackClass(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        projAnalysis(iTrack).momentScalingPowerSlope = mssSlope(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1);
        projAnalysis(iTrack).genDiffCoef = genDiffCoef(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1,:);
        projAnalysis(iTrack).momentScalingPower = scalingPower(...
            compTrackStartRow(iTrack):compTrackStartRow(iTrack)+...
            numSegments(iTrack)-1,:);
    end

end

%% ~~~ the end ~~~
