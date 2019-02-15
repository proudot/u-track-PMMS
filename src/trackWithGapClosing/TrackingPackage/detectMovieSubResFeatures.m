function detectMovieSubResFeatures(movieData,varargin)
% detectMovieSubResFeatures detect sub-resolution objects in a movie
%
% detectMovieSubResFeatures 
%
% SYNOPSIS detectMovieSubResFeatures(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   

% Sebastien Besson, Oct 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('SubResolutionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SubResolutionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
subResDetProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(subResDetProc,paramsIn);

%% --------------- Initialization ---------------%%

nChan=numel(movieData.channels_);
% Set up the input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
subResDetProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);
saveResults(nChan,1)=struct();
for i = p.ChannelIndex;    
    saveResults(i).dir = p.OutputDirectory ;
    saveResults(i).filename = ['Channel_' num2str(i) '_detection_result.mat'];
    %Create string for current directory
    outFilePaths{1,i} = [saveResults(i).dir filesep saveResults(i).filename ];
end
mkClrDir(p.OutputDirectory);
subResDetProc.setOutFilePaths(outFilePaths);

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting diffraction-limited objects...')

for i = p.ChannelIndex
    disp(['Please wait, detecting objects for channel ' num2str(i)])
    disp(inFilePaths{1,i});
    disp('Results will be saved under:')
    disp(outFilePaths{1,i});
    
    % Retrieve information about the images
    if movieData.isOmero() || movieData.isBF()
        movieParam.channel = movieData.channels_(i);
    else
        [~, base, digits4Enum, ~] = getFilenameBody(movieData.getImageFileNames{i}{1});
        digits4Enum = length(digits4Enum);
        
        movieParam.imageDir = [inFilePaths{1,i} filesep];
        movieParam.filenameBase = base;
        movieParam.digits4Enum = digits4Enum;
    end    
    movieParam.firstImageNum=p.firstImageNum;
    movieParam.lastImageNum=p.lastImageNum;

    % Call stand-alone subresolution detection function
    movieInfo = detectSubResFeatures2D_StandAlone(movieParam, p.detectionParam, saveResults(i));
end

disp('Finished detecting diffraction-limited objects...')

