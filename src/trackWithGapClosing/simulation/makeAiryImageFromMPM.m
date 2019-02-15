function [imageStack]=makeAiryImageFromMPM(trackInfo,bgav,bgnoise,sigma,imsize,rad,saveVar,saveFolder,saveDir)
% makeAiryImageFromMPM makes an image from the point distribution specified
% in the mpm file, using specified values for amplitudes and noise
%
% SYNOPSIS:
% [imageStack]=makeAiryImageFromMPM(mpm,amps,bg,bgnoise,sigma,imsize,rad)
%
% INPUT     :   trackInfo = MPM-like file with location and intensities of 
%                           objects in successive columns:
%                           x1,y1,i1,x2,y2,i2,...
%               bgav      = average background level (can be either single 
%                           value or, for background inequality, an image 
%                           of the same size as imsize). In counts, e.g.
%                           assuming a 16-bit camera.
%               bgnoise   = std of backgound noise. In counts, e.g.
%                           assuming a 16-bit camera.
%               sigma     = width of point-spread function (in pixels),
%                           where sigma = (1/3)*(radius of Airy disc)
%               imsize    = size of image [sx,sy]
%               rad       = radius used for the generation of the Airy disc
%                           for each object, in increments of sigma 
%                           ( should ideally be >=3 )
%               saveVar   = variable that indicates whether or not tif files
%                           are saved to file (1/0)
%               saveFolder  = (optional) folder name
%               saveDir   = Top directory to save images. The double
%                           "saveFolder" and "saveDir" are for historical
%                           reasons, for backward compatibility. In
%                           reality, saveFolder can be skipped and saveDir
%                           is enough. -KJ
%                           
%
% OUTPUT    :   imageStack
%
% REMARKS
%
% created with MATLAB ver.: 7.1.0.246 (R14) Service Pack 3 on Windows_NT
%
% created by: dloerke
% DATE: 28-Feb-2006

%% initialize parameters
ordir = cd;

[nx,ny]=size(trackInfo);
nframes = round(ny/3);

xs = imsize(1);
ys = imsize(2);

specImage = zeros(xs,ys);

% path to save images
saveYN = 0;
if nargin>6
    saveYN = saveVar;
    if saveYN == 1
        if nargin < 9 || isempty(saveDir)
            saveDir = uigetdir(ordir,'specify directory for saving simulation images');
        end
    end
end



%% generate little mask for individual airys
%generate list of pixels - circle of radius ns*sigma around each speckle,
%the value of ns is max(3,rad)
if (nargin>5)
    % if a value for rad is entered, use it, but it can't be less than the
    % Airy disc radius, which is specified by 3*sigma
    ns = max(rad,3);
else
    % default value is 6
    ns = 6;
end
% radius of cell mask in pixels
radius = ceil(ns*sigma);
% length of window
len = 2*radius+1;

% positions of mini-window in grid form
[miniImX, miniImY] = ndgrid(-radius:radius,-radius:radius);
miniDist = sqrt(miniImX.^2 + miniImY.^2);
% extract from the square grid positions those below the specified radius
% (convert to circular area)
[XinPix, YinPix] = find( miniDist<=radius );
XinPix = XinPix-(radius+1);
YinPix = YinPix-(radius+1);
% now XinPix and YinPix contain the positions of the pixels in the circular
% area around the origin with radius radius



%% make images

hr0 = waitbar(0,'image series','Position',[200 100 350 56.25]);

% loop over all frames
for t=1:nframes
    
    if saveYN==1
        cd(saveDir)
        
        % if a specific folder name is specified, go to this folder or
        % create it
        
        if nargin > 7 && ~isempty(saveFolder)
            % if folder of specified name already exists, go there
            if exist(saveFolder)==7
                cd(saveFolder);
            % else make it first and then go there
            else
                mkdir(saveFolder);
                cd(saveFolder);
            end
        end

            
        if t<10
            saveName = ['simulImage00',num2str(t),'.tif'];
        elseif t<100
            saveName = ['simulImage0',num2str(t),'.tif'];
        else
            saveName = ['simulImage',num2str(t),'.tif'];
        end
    end
    
    
    % object positions in current frames
    current_xy = trackInfo(:,3*t-2:3*t-1);
    
    current_pos = find((current_xy(:,1)>0) & (current_xy(:,2)>0));
    
    %instead of using grid of specified size, use input xy positions
    xyvecCon = [current_xy(current_pos,1) current_xy(current_pos,2) ];
    [numObjCon,nn]=size(xyvecCon);
    msposCon=xyvecCon;

    % amplitudes of objects are specified in the amps matrix 
    current_amps = trackInfo(:,3*t);
    ampvecCon = current_amps(current_pos);
    msposCon(:,3)=ampvecCon;
    
    % msposCon now contains positions and intensities of all objects in
    % this frame
    
    %======================================================================
    % Determine the positions of all the adjacent pixels relevant for each
    % object up until the specified distance
    %======================================================================
    %hr1 = waitbar(0,'generating sparse mask','Position',[345 250 350 56.25]);
    % loop over all objects
    for n=1:numObjCon
        % mspos are successive object positions; add these to the object
        % mask positions determined previously
        CurrCoorXCon = round(msposCon(n,1))+XinPix;
        CurrCoorYCon = round(msposCon(n,2))+YinPix;
        % bad posititons = those outside the image
        goodPos = find( (CurrCoorXCon>=1) & (CurrCoorYCon>=1) & (CurrCoorXCon<=xs) & (CurrCoorYCon<=ys) );
        badPos = find( (CurrCoorXCon<1) | (CurrCoorYCon<1) | (CurrCoorXCon>xs) | (CurrCoorYCon>ys) );
        subIndGoodPosCon = sub2ind([xs ys],CurrCoorXCon(goodPos),CurrCoorYCon(goodPos));
    
        CurrCoorXCon(badPos)=nan;
        CurrCoorYCon(badPos)=nan;
        
        % store the acceptable positions associated with this object
        singleObjectPositionsCon{n} = [CurrCoorXCon(goodPos),CurrCoorYCon(goodPos),sub2ind([xs ys],CurrCoorXCon(goodPos),CurrCoorYCon(goodPos))];
        
        %waitbar( (n / numObjCon),hr1 );
        
    end % of for n
    %close(hr1);

    % initialize: exact position of the current objects
    x0 = msposCon(:,1);
    y0 = msposCon(:,2);
    
    CurrFitimage = zeros(xs,ys);
    
    
    %======================================================================
    % Build up the image using the vicinity of each objects as determined
    % earlier
    %======================================================================
    
    %hr2 = waitbar(0,'building images','Position',[345 250 350 56.25]);
    % loop over all objects
    for z=1:numObjCon
        
        % xy positions of vicinity pixels, as well as index positions
        xyloc = singleObjectPositionsCon{z}(:,1:2);
        iloc = singleObjectPositionsCon{z}(:,3);
        
        % distances of vicinity pixels from current object
        zz=sqrt( (xyloc(:,1)-x0(z)).^2 + (xyloc(:,2)-y0(z)).^2 );
        % normalized psf intensities for all pixels
        M = psf2D_x0y0D(zz,3*sigma);
        % scale intensities with specified amplitude
        FittedCurveSin = ampvecCon(z) * M;
        % add current object pixels to current image
        CurrFitimage(iloc) = CurrFitimage(iloc)+ FittedCurveSin;
        
        %waitbar( (z / numObjCon),hr2 );
    end % of for z
    %close(hr2);
    
    % now add to the image the noisy background  
    specImage = uint16(round(bgav + CurrFitimage + bgnoise * randn(xs,ys)));
    
    
    if saveYN==1, imwrite(specImage,saveName,'tif'); end
    
    imageStack(:,:,t) = specImage;
    
    % display results
    if mod(t,10)==0
        imshow(specImage,[]); pause(0.01); 
    end
    
    waitbar( (t / nframes),hr0 );
end % of for t

close (hr0)

cd(ordir);

end % of function



%% ========================================================================
%  subfunctions
%  ========================================================================



function M = psf2D_x0y0D(dist,arad)

%PSF2D generates a 2D point spread function model (Airy disc)
%
% SYNOPSIS M = psf2D_x0y0(dist,arad)
%
% INPUT    dist    : distances from zero 
%          arad    : radius of airy disc = 3*mssigma
%
% OUTPUT   M       : filter mask (odd dimensions) representing the Airy disk
%                    the filter is normalized to sum(M(:)) = 1 
%
% NOTE     the function calculates the Airy disk radius according to 
%          R0 = 0.61 * lambda / NA
%          
% REMARK - M is not normalized here!
%
% Alexandre Matov, January 7th, 2003

% R0 = 0.61 * lambda / NA;
% R1 = 7.02 / 3.83 * R0;    % position of the second root of the Bessel function 

% dim = 2*ceil(R1/pixSize)+1;
% arad = 0.61*lambda/(NA*pixelsize)

% [x,y] = meshgrid(-ceil(R1/pixSize):ceil(R1/pixSize));
dist(dist==0)=0.0001;

% shift center by a 1/1000 pixel to avoid division-by-0
                               % in the next step
% ds = d*pixSize/R0*3.83;
ds = dist/arad*3.83;
psfs = (besselj(1,ds)./ds);
psf = psfs.*conj(psfs);

norms = (besselj(1,0.0001)./0.0001);
norm = norms.*conj(norms);
% normalization to value at zero;
M = psf / norm;                                 
%M = M/max(M(:));

end % of subfunction