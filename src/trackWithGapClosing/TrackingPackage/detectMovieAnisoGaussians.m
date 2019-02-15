function detectMovieAnisoGaussians(movieData,varargin)
% detectMovieAnisoGaussians detect objects by fitting anisotropic Gaussians
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

% Sebastien Besson, May 2012

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('AnisoGaussianDetectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(AnisoGaussianDetectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
detProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(detProc,paramsIn);

%% --------------- Initialization ---------------%%

nChan=numel(movieData.channels_);
% Set up the input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
detProc.setInFilePaths(inFilePaths);
    
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
detProc.setOutFilePaths(outFilePaths);

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting anistropic Gaussians...')

roiMask = movieData.getROIMask;

for i = p.ChannelIndex
    disp(['Please wait, detecting objects for channel ' num2str(i)])
    disp(inFilePaths{1,i});
    disp('Results will be saved under:')
    disp(outFilePaths{1,i});
    
    movieInfo(movieData.nFrames_,1)=struct('xCoord',[],'yCoord',[],...
        'amp',[],'sigmaX',[],'sigmaY',[],'theta',[],'bkg',[]);
    progressText(0,'Detecting anisotropic Gaussians');
    for j=1:movieData.nFrames_
        I=double(movieData.channels_(i).loadImage(j));
        if ~isempty(p.MaskProcessIndex)
            maskProc = movieData.processes_{p.MaskProcessIndex};
            mask = maskProc.loadChannelOutput(p.MaskChannelIndex,j);
        else
            mask = true(movieData.imSize_);
        end
        % Call stand-alone subresolution detection function
        movieInfo(j) = cometDetection(I, mask & roiMask(:,:,j), p.psfSigma,...
            'mode',p.mode,'kSigma',p.kSigma,'alpha',p.alpha,'minDist',p.minDist);
        progressText(j/movieData.nFrames_,'Detecting anisotropic Gaussians');
    end
    save(outFilePaths{1,i} ,'movieInfo');
end

disp('Finished detecting objects...')

