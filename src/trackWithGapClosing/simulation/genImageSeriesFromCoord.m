function genImageSeriesFromCoord(featCoord,imageParam,saveParam)

%get information for saving movie
saveDir = saveParam.saveDir;
movieNameBase = saveParam.movieNameBase;

%define image parameters
psfSigma = imageParam.psfSigma;
imSize = imageParam.imSize;
signalAboveBG = imageParam.signalAboveBG;
SNR = imageParam.SNR;

%assign background mean and derive background std from image parameters
bgMean = 30000;
bgStd = signalAboveBG ./ SNR;

%get number of features and number of frames
[numFeat,~,numFrames] = size(featCoord);

%separate feature coordinates into integer + fraction
featCoordFloor = floor(featCoord);
featCoordFrac = featCoord - featCoordFloor;

%make a PSF template to get PSF range
template = signalAboveBG * GaussMask2D(psfSigma,6*ceil(psfSigma)-1,[0 0]);
psfRange = size(template,1);
psfRange = floor(psfRange/2);

%make "movie"...

%go over all frames
for iFrame = 1 : numFrames
    
    %make empty image
    imageFeatures = zeros(imSize+2*psfRange,imSize+2*psfRange);
    
    %go over all features
    for iFeat = 1 : numFeat
        if ~isnan(featCoordFrac(iFeat,:,iFrame))
            
            %get PSF template for this feature
            template = signalAboveBG * GaussMask2D(psfSigma,6*ceil(psfSigma)-1,featCoordFrac(iFeat,:,iFrame));
            
            x0 = featCoordFloor(iFeat,1,iFrame) + psfRange;
            y0 = featCoordFloor(iFeat,2,iFrame) + psfRange;
            
            xmin = x0 - psfRange;
            xmax = x0 + psfRange;
            
            ymin = y0 - psfRange;
            ymax = y0 + psfRange;
            
            %add feature PSF to image
            imageFeatures(xmin:xmax,ymin:ymax) = imageFeatures(xmin:xmax,ymin:ymax) + template;
            
        end
    end
    imageFeatures = imageFeatures(psfRange+1:end-psfRange,psfRange+1:end-psfRange);
    
    %pick random numbers for background
    randNumBG = randn(imSize,imSize);
    
    %construct images with background
    imageWithNoise = imageFeatures + bgMean + randNumBG*bgStd;
    
    %convert to uint16
    imageWithNoise = uint16(imageWithNoise);
    
    %save images in directory
    if iFrame < 10
        filename = [saveDir movieNameBase '00' num2str(iFrame) '.tif'];
        imwrite(imageWithNoise,filename,'tif','compression','none');
    elseif iFrame < 100
        filename = [saveDir movieNameBase '0' num2str(iFrame) '.tif'];
        imwrite(imageWithNoise,filename,'tif','compression','none');
    else
        filename = [saveDir movieNameBase num2str(iFrame) '.tif'];
        imwrite(imageWithNoise,filename,'tif','compression','none');
    end

end

