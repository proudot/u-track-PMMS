
%% general gap closing parameters
gapCloseParam.timeWindow = 1; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatStationaryLink';

%parameters
parameters.searchRadius = 10;

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatStationaryCloseGaps';

%parameters
parameters.searchRadius = 10; 
parameters.gapPenalty = 2;

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions = [];

%% additional input

%saveResults
saveResults.dir = '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110825_FixedCellTIRF_CD36/10000/Field_5/Analysis/'; %directory where to save input and output
saveResults.filename = 'tracksDetectionAll2.mat'; %name of file where input and output are saved
% saveResults = 0; %don't save results

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% tracking function call

[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

% for startFrame = 1 : 400 : 48000
%     endFrame = startFrame + 399;
%     saveResults.filename = ['tracks2Detection1_Frames' sprintf('%05i',startFrame) 'to' sprintf('%05i',endFrame) '.mat'];
%     disp(startFrame)
%     [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(...
%         movieInfo(startFrame:endFrame),costMatrices,gapCloseParam,...
%         kalmanFunctions,probDim,saveResults,verbose);
% end

%% ~~~ the end ~~~