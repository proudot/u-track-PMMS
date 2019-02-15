% featuresInfo = cometDetection(img, mask, psfSigma, mode)
%
% Inputs :      img : input image
%              mask : cell mask
%             psfSigma : standard deviation of the Gaussian PSF
%            {mode} : parameters to estimate, default 'xyArtc'
%           {alpha} : alpha-values for statistical test
%           {kSigma} : alpha-values for statistical test
%           {minDist} : minimum distance betwen detected features
%           {filterSigma} : sigma for the steerable filter
%
% Outputs:  featuresInfo : output structure with anisotropic Gaussian
%                          parameters, standard deviations (compatible with
%                          Khuloud's tracker.
%
% Sylvain Berlemont, April 2011

function featuresInfo = cometDetection(img, mask, psfSigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('mask', @islogical);
ip.addRequired('psfSigma', @isscalar);
ip.addParamValue('mode', 'xyArtc', @ischar);
ip.addParamValue('alpha', 0.05, @isscalar);
ip.addParamValue('kSigma', 4, @isscalar);
ip.addParamValue('minDist', .25, @isscalar);
ip.addParamValue('filterSigma',psfSigma*sqrt(2), @isscalar);
ip.addParamValue('offset','centered', @ischar);
ip.addParamValue('threshold','otsuRosin', @ischar);


ip.parse(img, mask, psfSigma, varargin{:});
mode = ip.Results.mode;
alpha = ip.Results.alpha;
kSigma = ip.Results.kSigma;
minDist = ip.Results.minDist;
filterSigma =ip.Results.filterSigma;

img = double(img);
[nrows ncols] = size(img);

% Filter image with laplacian
bandPassIso = filterLoG(img,psfSigma);
bandPassIso(bandPassIso < 0) = 0;
bandPassIso(~mask) = 0;

% Filter image with steerable filter
[R,T] = steerableDetector(img,2,filterSigma);
 
% Compute the local maxima of the bandpass filtered images
locMaxIso = locmax2d(R, [5 5]);

bw=[];
switch ip.Results.threshold
    case 'otsuRosin'
        bw = blobSegmentThreshold(bandPassIso,0,0,mask);
    case 'pointSource'
        [bw,imgLM] = pointSourceStochasticFiltering2D(img, filterSigma);
    otherwise 
        error('Unknown thresholding method.');
end

% figure()
% subplot(1,2,1); imshow(blobSegmentThreshold(bandPassIso,0,0,mask));
% subplot(1,2,2); imshow(pointSourceStochasticFiltering2D(img, filterSigma));

labels=bwlabel(bw);

locMaxIso(~bw) = 0;

indMax = find(locMaxIso);
[y x] = ind2sub(size(img), indMax);

P = zeros(size(y, 1), 7);
P(:,1) = x;
P(:,2) = y;
P(:,3) = img(indMax);
P(:,4) = 2*psfSigma;       % sigmaX  % !!!! suppress 2x (repeated below)
P(:,5) = psfSigma;         % sigmaY
P(:,6) = T(indMax)+pi/2;

% % Subresolution detection
% hside = ceil(kSigma * psfSigma);
% npx = (2 * hside + 1)^2;
% xmin = x - hside;
% xmax = x + hside;
% ymin = y - hside;
% ymax = y + hside;

PP =num2cell(P,1);

[xRange,yRange,nzIdx] = arrayfun(@(x0,y0,sigmaX,sigmaY,theta)...
    anisoGaussian2DSupport(x0,y0,sigmaX,sigmaY,theta,kSigma,[ncols nrows]),...
    PP{[1 2 4 5 6]},'UniformOutput',false);
% hside = ceil(kSigma * psfSigma);
% npx = (2 * hside + 1)^2;
xmin = cellfun(@min,xRange);
xmax = cellfun(@max,xRange);
ymin = cellfun(@min,yRange);
ymax = cellfun(@max,yRange);
npx = cellfun(@numel,nzIdx);

isValid = find(xmin >= 1 & xmax <= ncols & ymin >= 1 & ymax <= nrows);

xmin = xmin(isValid);
xmax = xmax(isValid);
ymin = ymin(isValid);
ymax = ymax(isValid);
P = P(isValid,:);

stdP = zeros(size(P));
stdR = zeros(size(P,1),1);

kLevel = norminv(1 - alpha / 2.0, 0, 1); % ~2 std above background

success = false(numel(xmin),1);

for iFeature = 1:numel(xmin)
    mask =labels(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    labelMask=mask;
    % Set to zero the pixel that belong to the central object
    centralLabel=labels(yRange{iFeature}(fix(end/2)),xRange{iFeature}(fix(end/2)));
    if(centralLabel~=0)
        mask(mask==centralLabel) = 0; 
    else
        mask=0;
    end
        
    %%
    anisoMask = false(length(yRange{iFeature}),length(xRange{iFeature}));
    anisoMask(nzIdx{iFeature})=true;
    %%
    anisoMask(mask~=0)=false; % pixel of the other label are suppressed
    %%
    npx(iFeature)=nnz(anisoMask);
    crop = img(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    cropNoMask = crop;
    crop(~anisoMask)=NaN;

    
    P(iFeature,7) = min(crop(:)); % background
    P(iFeature,3) = P(iFeature,3) - P(iFeature,7); % amplitude above background
        
    [params, stdParams, ~, res] = fitAnisoGaussian2D(crop, ...
        [0, 0, P(iFeature,3),3*  P(iFeature,4), P(iFeature,5), ...
        P(iFeature,6), P(iFeature,7)], mode);
        
    % TEST: position must remain in a confined area
    px = floor(floor(size(crop)/2)+1+params(1:2));
    isValid = all(px>=1) & all(px<=size(crop));
    if isValid
        isValid = anisoMask(px(1),px(2));
    end
    
    % TEST: sigmaX > 1
    isValid = isValid & params(4) > 1;
    
%     TEST: goodness-of-fit
%     stdRes = std(res.data(~isnan(res.data)));
%     [~, pval] = kstest(res.data(~isnan(res.data)) ./ stdRes, [], alpha);
    isValid = isValid & res.pval > alpha;

    % TEST: amplitude
    SE_psfSigma_r = (res.std / sqrt(2*(npx(iFeature)-1))) * kLevel;
    psfSigma_A = stdParams(3);
    A_est = params(3);
    df2 = (npx(iFeature) - 1) * (psfSigma_A.^2 + SE_psfSigma_r.^2).^2 ./ (psfSigma_A.^4 + SE_psfSigma_r.^4);
    scomb = sqrt((psfSigma_A.^2 + SE_psfSigma_r.^2) / npx(iFeature));
    T = (A_est - res.std * kLevel) ./ scomb;    
    isValid = isValid & (1 - tcdf(T, df2)) < alpha;
  
    % TEST: extreme value of 
    isValid = isValid & params(4) < 10 * psfSigma;
    
    success(iFeature) = isValid;
    
    LMX=P(iFeature,1); 
    LMY=P(iFeature,2); 
    
    P(iFeature,1) = P(iFeature,1) + params(1);
    P(iFeature,2) = P(iFeature,2) + params(2);
    P(iFeature,3) = params(3);
    P(iFeature,4) = params(4);
    P(iFeature,5) = params(5);
    P(iFeature,6) = params(6);
    P(iFeature,7) = params(7);
   
    cropMasked=crop;cropMasked(~anisoMask)=0;
    cropMasked(isnan(cropMasked))=0;
    switch ip.Results.offset
        case 'tipCentroid'
            %%
            maskOne=labelMask;maskOne(maskOne>1)=1;
            prop = regionprops(maskOne, (cropNoMask), {'Centroid','WeightedCentroid'});
            wCentroid=prop(1).WeightedCentroid-[(size(crop,2)/2) (size(crop,1)/2)];
            P(iFeature,1) = wCentroid(1) + LMX;
            P(iFeature,2) = wCentroid(2) + LMY;
            %%
            showDebug=0;
            if(showDebug)
                subplot(1,4,1)
                imshow(mat2gray(crop))
                hold on
                scatter((size(crop,2)/2),(size(crop,1)/2),'b+');
                scatter((size(crop,2)/2)+params(1),(size(crop,1)/2)+params(2),'ro');
                scatter((size(crop,2)/2)+wCentroid(1),(size(crop,1)/2)+ wCentroid(2),'go');
                hold off
                subplot(1,4,2)
                imshow(mat2gray(labelMask))
                subplot(1,4,3)
                imshow(mat2gray(cropMasked.^2))
                
                
            end
            
        case 'tipFit'
            %% Integrate value along the seam
            theta=params(6);
            cropLabelMask=crop;cropLabelMask(labelMask==0)=0;
            rotcrop=imrotate(cropLabelMask,180*(theta/pi));
            rotcrop(isnan(rotcrop))=0;
            prof=sum(rotcrop)./sum(rotcrop>0);            
            prof(isnan(prof))=0;
          % prof=smooth(prof,3);
            [m,i]=max(prof);
            posCentered=i-length(prof)/2;
            
            R = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];
            rotCoord=R*[params(2); params(1)];
            rotTipCoord=[posCentered; rotCoord(2)];
            RBack = inv(R);
            tipCoordCrop=RBack*rotTipCoord;
            rotCoordInv=RBack*rotCoord;
            
            
            showDebug=0;
            if(showDebug)
                subplot(1,5,1)
                imshow(mat2gray(crop))
                hold on
                scatter((size(crop,2)/2),(size(crop,1)/2),'b+');
                scatter((size(crop,2)/2)+params(1),(size(crop,1)/2)+params(2),'ro');
                scatter((size(crop,2)/2)+tipCoordCrop(1),(size(crop,1)/2)+tipCoordCrop(2),'go');
                hold off
                subplot(1,5,2)
                imshow(mat2gray(labelMask))
                 hold on
                scatter(size(labelMask,2)/2,size(labelMask,1)/2,'g+')
                hold off 
                subplot(1,5,3)
                imshow(mat2gray(cropNoMask))
                                hold on
                scatter(size(cropNoMask,2)/2,size(cropNoMask,1)/2,'g+')
                hold off 
                subplot(1,5,4)
                imshow(mat2gray(rotcrop))
                                hold on
                scatter((size(rotcrop,2)/2),(size(rotcrop,1)/2),'b+');  
                scatter((size(rotcrop,2)/2)+rotCoord(1),(size(rotcrop,1)/2)+rotCoord(2),'ro');
                scatter((size(rotcrop,2)/2)+rotTipCoord(1),(size(rotcrop,1)/2)+rotTipCoord(2),'go');
                hold off
                subplot(1,5,5)
                plot(prof)
            end
            
             P(iFeature,1) = LMX+tipCoordCrop(1);
             P(iFeature,2) = LMY+tipCoordCrop(2);     
        case 'center'     
            showDebug=0;
            if(showDebug)
                subplot(1,4,1)
                imshow(mat2gray(crop))
                hold on
                scatter((size(crop,2)/2),(size(crop,1)/2),'b+');
                scatter((size(crop,2)/2)+params(1),(size(crop,1)/2)+params(2),'ro');
                hold off
                subplot(1,4,2)
                imshow(mat2gray(labelMask))
                subplot(1,4,3)
                imshow(mat2gray(cropMasked.^2))               
            end
        otherwise
            error('Unknown method');
    end
    
    stdP(iFeature,1) = stdParams(1);
    stdP(iFeature,2) = stdParams(2);
    stdP(iFeature,3) = stdParams(3);
%     stdP(iFeature,4) = stdParams(4);
    stdP(iFeature,6) = stdParams(4);
    stdP(iFeature,7) = stdParams(5);
    
    stdR(iFeature) = res.std;
end

P = P(success,:);
stdP = stdP(success,:);

% Remove any detection which has been localised at the same position
isValid = true(size(P,1),1);
idxKD = KDTreeBallQuery(P(:,1:2), P(:,1:2), repmat(minDist, size(P,1), 1));
idxKD = idxKD(cellfun(@(x) length(x)>1, idxKD));
    
for k = 1:length(idxKD);
    stdRes = stdR(idxKD{k});
    isValid(idxKD{k}(stdRes ~= min(stdRes))) = false;
end

P = P(isValid,:);
stdP = stdP(isValid,:);

featuresInfo.xCoord = [P(:,1), stdP(:,1)];
featuresInfo.yCoord = [P(:,2), stdP(:,2)];
featuresInfo.amp = [P(:,3), stdP(:,3)];
featuresInfo.sigmaX = [P(:,4), stdP(:,4)];
featuresInfo.sigmaY = [P(:,5), stdP(:,5)];
featuresInfo.theta = [P(:,6), stdP(:,6)];
featuresInfo.bkg = [P(:,7), stdP(:,7)];
