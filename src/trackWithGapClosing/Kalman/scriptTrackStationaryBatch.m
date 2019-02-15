%% define batch job locations

%movieInfo locations
movieInfoDir = {...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_7/Field1/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_7/Field2/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_6/Field1/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_6/Field2/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_5/Field1/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_5/Field2/Analysis/',...
    };

%movieInfo file names
movieInfoFile = {...
    'detectionAll2.mat',...
    'detectionAll2.mat',...
    'detectionAll2.mat',...
    'detectionAll2.mat',...
    'detectionAll2.mat',...
    'detectionAll2.mat',...
    };

%directory for saving results
saveResDir = {...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_7/Field1/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_7/Field2/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_6/Field1/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_6/Field2/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_5/Field1/Analysis/',...
    '/home/kj35/files/LCCB/receptors/endoCells/AnindyaChanda/20110822_Monodispersion/10_5/Field2/Analysis/',...
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
for iMovie = 1 : numMovies1
    
    try
        
        %display movie number
        disp(['Movie ' num2str(iMovie) ' / ' num2str(numMovies1) ' ...'])
        
        %get movieInfo
        load(fullfile(movieInfoDir{iMovie},movieInfoFile{iMovie}));
        
        %% general gap closing parameters
        gapCloseParam.timeWindow = 200; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
        gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
        gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.
        
        %optional input:
        gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
        
        %% cost matrix for frame-to-frame linking
        
        %function name
        costMatrices(1).funcName = 'costMatStationaryLink';
        
        %parameters
        parameters.searchRadius = 1.5;
        
        costMatrices(1).parameters = parameters;
        clear parameters
        
        %% cost matrix for gap closing
        
        %function name
        costMatrices(2).funcName = 'costMatStationaryCloseGaps';
        
        %parameters
        parameters.searchRadius = 1.5;
        parameters.gapPenalty = 2;
        
        costMatrices(2).parameters = parameters;
        clear parameters
        
        %% Kalman filter function names
        
        kalmanFunctions = [];
        
        %% additional input
        
        %saveResults
        saveResults.dir = saveResDir{iMovie}; %directory where to save input and output
        saveResults.filename = 'tracksDetectionAll2.mat'; %name of file where input and output are saved
        % saveResults = 0; %don't save results
        
        %verbose
        verbose = 1;
        
        %problem dimension
        probDim = 2;
        
        %% tracking function call
        
        [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
            costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);
        
        %         for startFrame = 1 : 400 : 6800
        %             endFrame = startFrame + 399;
        %             saveResults.filename = ['tracks8Detection2_Frames' sprintf('%04i',startFrame) 'to' sprintf('%04i',endFrame) '.mat'];
        %             disp(startFrame)
        %             [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(...
        %                 movieInfo(startFrame:endFrame),costMatrices,gapCloseParam,...
        %                 kalmanFunctions,probDim,saveResults,verbose);
        %         end
        
    catch %#ok<CTCH>
        disp(['Movie ' num2str(iMovie) ' failed!']);
    end
    
end
