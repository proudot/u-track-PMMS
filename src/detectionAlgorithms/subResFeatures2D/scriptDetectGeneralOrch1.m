
%% movie information
movieParam.imageDir = '/files/LCCB/receptors/Galbraiths/data/beta3andCellEdge/130503_Cs2C1_Beta3/imagesBeta3/'; %directory where images are
movieParam.filenameBase = '130503_Cs2C1_mEos2Beta3_'; %image file name base
movieParam.firstImageNum = 2; %number of first image in movie
movieParam.lastImageNum = 14400; %number of last image in movie
movieParam.digits4Enum = 5; %number of digits used for frame enumeration (1-4).

%% detection parameters
detectionParam.psfSigma = 1.2; %point spread function sigma (in pixels)
detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.05,'alphaD',0.05,'alphaF',0); %alpha-values for detection statistical tests
detectionParam.visual = 0; %1 to see image with detected features, 0 otherwise
detectionParam.doMMF = 1; %1 if mixture-model fitting, 0 otherwise
detectionParam.bitDepth = 16; %Camera bit depth
detectionParam.alphaLocMax = 0.1; %alpha-value for initial detection of local maxima
detectionParam.numSigmaIter = 0; %maximum number of iterations for PSF sigma estimation
detectionParam.integWindow = 0; %number of frames before and after a frame for time integration

detectionParam.calcMethod = 'g';

%absolute background info and parameters...
background.imageDir = '/files/LCCB/receptors/Galbraiths/data/beta3andCellEdge/130503_Cs2C1_Beta3/bgBeta3/';
background.filenameBase = 'crop_130503_Cs2C1_mEos2Beta3_';
background.alphaLocMaxAbs = 0.01;
detectionParam.background = background;

%% additional input

%save results
saveResults.dir = '/files/LCCB/receptors/Galbraiths/data/beta3andCellEdge/130503_Cs2C1_Beta3/'; %directory where to save input and output
saveResults.filename = 'detectionAll1.mat'; %name of file where input and output are saved
% saveResults = 0;

%verbose state
verbose = 0;

%% run the detection function
[movieInfo,exceptions,localMaxima,background,psfSigma] = ...
    detectSubResFeatures2D_StandAlone(movieParam,detectionParam,saveResults,verbose);

%% exit

exit
