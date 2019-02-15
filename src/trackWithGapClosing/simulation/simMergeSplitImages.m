% Initialize image dimensions
nx = 256;
ny = 256;
sizeT = 50;

dx = 2;
dy = 3;
x0 = nx/4;

%% Create 3 merging and splitting tracks

% Master track
x1 = x0:dx:x0+(sizeT - 1)*dx;
y1 = ny/2 * ones(1, sizeT);

% Second track (merging at sizeT/5 and splitting at 2sizeT/5)
x2 = zeros(size(x1));
y2 = zeros(size(x1));
midrange = sizeT/5 + 1 : 2 * sizeT/5;
x2(midrange) = x1(midrange);
y2(midrange) = y1(midrange);

firstrange = 1 : sizeT/5;
x2(firstrange) = x2(midrange(1));
y2(firstrange) = y2(midrange(1)) + dy * (length(firstrange)-1:-1:0);

endrange = 2 * sizeT/5+1 : sizeT;
x2(endrange) = x2(midrange(end));
y2(endrange) = y2(midrange(end)) - dy * (0:length(endrange)-1);

% Third track (merging at 3*sizeT/5 and splitting at 4*sizeT/5)
x3 = zeros(size(x1));
y3 = zeros(size(x1));
midrange = 3 * sizeT/5 + 1 : 4 * sizeT/5;
x3(midrange) = x1(midrange);
y3(midrange) = y1(midrange);

firstrange = 1 : 3 * sizeT/5;
x3(firstrange) = x3(midrange(1));
y3(firstrange) = y3(midrange(1)) +  dy * (length(firstrange)-1:-1:0);

endrange = 4 * sizeT/5+1 : sizeT;
x3(endrange) = x3(midrange(end));
y3(endrange) = y3(midrange(end)) - dy * (0:length(endrange)-1);

x = [x1; x2; x3];
y = [y1; y2; y3];
amp = 20 * ones(size(x));

%% Create synthetic image of merge and splitting

% Initialize array
I = 10 * rand(ny, nx, 1, 1, sizeT) + 10;

sigma=1.4;
% Create gap index for master track
gapInterval = sizeT/2 -1 : sizeT/2; 
validTracks = true(size(x));
validTracks(1, gapInterval) = false;

% Create diffraction limited spots
for t = 1 : sizeT
    I(:, :, 1, 1, t) = I(:, :, 1, 1, t) + simGaussianSpots(nx, ny, sigma,...
        'x', x(validTracks(:,t), t), 'y', y(validTracks(:,t), t),...
        'A', amp(validTracks(:,t), t), 'npoints', sum(validTracks(:,t)));        
end

% Save synthetic image as OME-TIFF
outFile = '/tmp/mergesplit.ome.tiff';
if exist(outFile, 'file'), delete(outFile); end
bfsave(uint8(I), outFile);

%% Analysis

% Initialize MovieData and set basic metadata
MD = MovieData.load(outFile);
MD.pixelSize_ = 67;
MD.getChannel(1).emissionWavelength_ = 500;
MD.numAperture_ = 1.4;
MD.sanityCheck();

% Set up tracking package
MD.addPackage(UTrackPackage(MD));
MD.getPackage(1).createDefaultProcess(1);
MD.getPackage(1).createDefaultProcess(2);
p = MD.getPackage(1).getProcess(2).getParameters();
p.gapCloseParam.mergeSplit = 1;
MD.getPackage(1).getProcess(2).setParameters(p);

% Run detection and tracking
MD.getPackage(1).getProcess(1).run();
MD.getPackage(1).getProcess(2).run();