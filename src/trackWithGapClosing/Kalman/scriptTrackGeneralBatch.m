%% define batch job locations

%movieInfo locations
movieInfoDir = {...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131127\analysisAlphaVY773A\',...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131202\analysisAlphaVY773A\',...
    };

%movieInfo file names
movieInfoFile = {...
    'detectionAll1.mat',...
    'detectionAll1.mat',...
    };

%directory for saving results
saveResDir = {...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131127\analysisAlphaVY773A\',...
    'C:\kjData\Galbraiths\data\alphaVY773AandCellEdge\131202\analysisAlphaVY773A\',...
    };

%calculate number of movies
numMovies1 = length(movieInfoDir);
numMovies2 = length(movieInfoFile);
numMovies3 = length(saveResDir);
if any([numMovies1-numMovies2 numMovies2-numMovies3 numMovies3-numMovies1]) ~= 0
    disp('Number of entries in movieInfoDir, movieInfoFile and saveResDir does not match!')
    return
end

%% track
% for iMovie = 1 : numMovies1
for iMovie = 2 : 2
    
    try
        
        %display movie number
        disp(['Movie ' num2str(iMovie) ' / ' num2str(numMovies1) ' ...'])
        
        %get movieInfo
        load(fullfile(movieInfoDir{iMovie},movieInfoFile{iMovie}));
        
        %% general gap closing parameters
        gapCloseParam.timeWindow = 7; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
        gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
        gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
        
        %optional input:
        gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
        
        %% cost matrix for frame-to-frame linking
        
        %function name
        costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';
        
        %parameters
        
        parameters.linearMotion = 0; %use linear motion Kalman filter.
        
        parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
        parameters.maxSearchRadius = 4.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
        parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
        
        parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
        parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
        % parameters.nnWindow = 10; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
        
        parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
        
        %optional input
        parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.
        
        costMatrices(1).parameters = parameters;
        clear parameters
        
        %% cost matrix for gap closing
        
        %function name
        costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';
        
        %parameters
        
        %needed all the time
        parameters.linearMotion = 0; %use linear motion Kalman filter.
        
        parameters.minSearchRadius = 2; %minimum allowed search radius.
        parameters.maxSearchRadius = 4.5; %maximum allowed search radius.
        parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
        
        parameters.brownScaling = [0.25 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
        parameters.timeReachConfB = gapCloseParam.timeWindow; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
        
        parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
        
        parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
        
        parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
        parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
        % parameters.nnWindow = 10; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
        
        parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
        
        parameters.linScaling = [0.25 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
        parameters.timeReachConfL = gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
        
        parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
        
        %optional; if not input, 1 will be used (i.e. no penalty)
        parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
        
        %optional; to calculate MS search radius
        %if not input, MS search radius will be the same as gap closing search radius
        parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.
        
        costMatrices(2).parameters = parameters;
        clear parameters
        
        %% Kalman filter function names
        
        kalmanFunctions.reserveMem  = 'kalmanResMemLM';
        kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
        kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
        kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
        
        %% additional input
        
        %saveResults
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        %         saveResults.filename = 'tracksFinalDetection3.mat'; %name of file where input and output are saved
        % saveResults = 0; %don't save results
        
        %verbose
        verbose = 1;
        
        %problem dimension
        probDim = 2;
        
        %% tracking function call
        
        %         [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
        %             costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
        
        maxI = length(movieInfo)/1200;
        %         for i = 1 : maxI
        for i = 2 : maxI
            disp(num2str(i))
            movieInfoTmp((i-1)*1200+1:i*1200) = movieInfo((i-1)*1200+1:i*1200);
            saveResults.filename = ['tracks1All_' sprintf('%02i',i) '.mat'];
            [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(...
                movieInfoTmp,costMatrices,gapCloseParam,...
                kalmanFunctions,probDim,saveResults,verbose);
        end
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
