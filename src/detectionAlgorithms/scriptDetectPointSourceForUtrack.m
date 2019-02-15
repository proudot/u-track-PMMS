
%% movie information
movieParam.imageDir = '../data/';
movieParam.filenameBase = 'cell042_'; %image file name base
movieParam.firstImageNum = 1; %number of first image in movie
movieParam.lastImageNum = 61; %number of last image in movie
movieParam.digits4Enum = 4; %number of digits used for frame enumeration (1-4).

%% detection parameters
pathsList=getImgPaths(movieParam);
estimated_sigma=getGaussianPSFsigmaFromData(pathsList(1:20))

sigma=estimated_sigma;
alpha=0.00001; %         'alpha' : alpha value used in the statistical tests. Default: 0.05.
verbose = 1;


%% additional input

%saveResults
saveResults.dir =  ['../res']; %directory where to save input and output
mkdir(saveResults.dir);
saveResults.filename = 'detectionPointSource.mat'; %name of file where input and output are saved


%% run the detection function
[movieInfo] = pointSourceDetect_StandAlone(movieParam,saveResults,verbose,sigma,8,'alpha',alpha);


%% display
startend=[movieParam.firstImageNum movieParam.lastImageNum];
saveMovie=1;
movieName='DetectionPointSource';
filterSigma=0;
showRaw=1;
intensityScale=0;
firstImageFile=strcat(movieParam.imageDir,movieParam.filenameBase,...
											sprintf('%04d ',movieParam.firstImageNum),'.tif');
dir2saveMovie=saveResults.dir;


overlayFeaturesMovie(movieInfo,[movieParam.firstImageNum movieParam.lastImageNum],saveMovie,movieName,...
										 filterSigma,showRaw,intensityScale, ...
                     firstImageFile,dir2saveMovieminLength);


