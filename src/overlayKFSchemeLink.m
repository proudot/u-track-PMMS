function overlayKFSchemeLink(movieInfo,tracksFeatIndxLink,KalmanInfoLink,costMatLinkParam,startend,dragtailLength,...
                             saveMovie,movieName,filterSigma,classifyGaps,highlightES,showRaw,...
                             imageRange,onlyTracks,classifyLft,diffAnalysisRes,intensityScale,...
                             colorTracks,firstImageFile,dir2saveMovie,varargin)
%OVERLAY 
%
%SYNPOSIS In this section we do not think in track but in link. For
%         each track (or segment) in tracksFeatIndxLink each segment is
%         consider and sorted in a list correponding to its associated KF
%         scheme.

%INPUT  tracksFinal   : Output of trackCloseGapsKalman.
%       startend      : Row vector indicating first and last frame to
%                       include in movie. Format: [startframe endframe].
%                       Optional. Default: [(first frame with tracks) (last frame with tracks)]
%       dragtailLength: Length of drag tail (in frames).
%                       Optional. Default: 10 frames.
%                       ** If dragtailLength = 0, then no dragtail.
%                       ** To show tracks from their beginning to their end,
%                       set dragtailLength to any value longer than the
%                       movie.
%                       ** To show tracks statically while features dance
%                       on them, use -1.
%                       ** To show tracks from their beginning to their
%                       end, and to retain tracks even after the particle
%                       disappears, use -2.
%       saveMovie     : 1 to save movie (as Quicktime), 0 otherwise.
%                       Optional. Default: 0.
%       movieName     : filename for saving movie.
%                       Optional. Default: TrackMovie (if saveMovie = 1).
%       filterSigma   : 0 to overlay on raw image, PSF sigma to overlay on
%                       image filtered with given filterSigma.
%                       Optional. Default: 0.
%       classifyGaps  : 1 to classify gaps as "good" and "bad", depending
%                       on their length relative to the legnths of the
%                       segments they connect, 0 otherwise.
%                       Optional. Default: 1.
%       highlightES   : 1 to highlight track ends and starts, 0 otherwise.
%                       Optional. Default: 1.
%       showRaw       : 1 to add raw movie to the left of the movie with
%                       tracks overlaid, 2 to add raw movie at the top of
%                       the movie with tracks overlaid, 0 otherwise.
%                       Optional. Default: 0.
%       imageRange    : Image region to make movie out of, in the form:
%                       [min pixel X, max pixel X; min pixel Y, max pixel Y].
%                       Optional. Default: Whole image.
%       onlyTracks    : 1 to show only tracks without any symbols showing
%                       detections, closed gaps, merges and splits; 0 to
%                       show symbols on top of tracks.
%                       Optional. Default: 0.
%       classifyLft   : 1 to classify objects based on (1) whether they
%                       exist throughout the whole movie, (2) whether they
%                       appear OR disappear, and (3) whether they appear
%                       AND disappear; 0 otherwise.
%                       Optional. Default: 1.
%       diffAnalysisRes:Diffusion analysis results (either output of
%                       trackDiffusionAnalysis1 or trackTransientDiffusionAnalysis2).
%                       Needed if tracks/track segments are to be
%                       colored based on their diffusion classification.
%                       With this option, classifyGaps, highlightES and
%                       classifyLft are force-set to zero, regardless of input.
%                       Optional. Default: None.
%       intensityScale: 0 to autoscale every image in the movie, 1
%                       to have a fixed scale using intensity mean and std,
%                       2 to have a fixed scale using minimum and maximum
%                       intensities.
%                       Optional. Default: 1.
%       colorTracks   : 1 to color tracks by rotating through 7 different
%                       colors, 0 otherwise. With this option,
%                       classifyGaps, highlightES and classifyLft are
%                       force-set to zero, regardless of input.
%                       Option ignored if diffAnalysisRes is supplied.
%                       Optional. Default: 0.
%       firstImageFile: Name of the first image file in the folder of
%                       images that should be overlaid. The file has to be
%                       the first image that has been analyzed even if not
%                       plotted. If file is not specified [], user will be
%                       prompted to select the first image.
%                       Optional. Default: [].
%       dir2saveMovie:  Directory where to save output movie.
%                       If not input, movie will be saved in directory where
%                       images are located.
%                       Optional. Default: [].
%       minLength     : Minimum length of tracks to be ploted.
%                       Optional. Default: 1.
%       plotFullScreen: 1 the figure will be sized to cover the whole
%                       screen. In this way the movie will be of highest
%                       possible quality. default is 0.
%       movieType     : 'mov' to make a Quicktime movie using MakeQTMovie,
%                       'avi' to make AVI movie using Matlab's movie2avi,
%                       'mp4_unix', 'avi_unix' to make an MP4 or AVI movie
%                       using ImageMagick and ffmpeg. These options work
%                       only under linux or mac.
%                       Optional. Default: 'mov'.
%
%OUTPUT the movie.
%
%REMARKS Color-coding:
%        ** Without diffusion classification, all tracks have a neutral
%        color, while objects are color coded in the following way:
%               * Detected object just after appearance: Green circle.
%               * Detected object just before disappearance: Yellow
%                 circle.
%               * Detected object in middle of trajectory that spans
%                 whole movie: White circle.
%               * Detected object in middle of trajectory that appears OR
%                 disappears within movie: Magenta circle.
%               * Detected object in middle of trajectory that appears AND
%                 disappears within movie: Red circle.
%               * Gap that is short than both segments it connects: Cyan
%                 star.
%               * Gap that is longer than at least one ofthe segments it
%                 connects: Blue star.
%               * Object before and after splitting: Green diamond.
%               * OBject before and after merging: Yellow diamond.
%           When classifyGaps = 0, all gaps are cyan.
%           When highlighES = 0, no green and yellow circles.
%           When classifyLft = 0, all objets in middle of trajectory are white.
%
%       ** With diffusion classification, all objects and gaps have neutral
%       color (merges and splits are diamonds), while tracks and track
%       segments are color-coded in the following way:
%               * Type 1: Linear + 1D confined diffusion: Orange.
%               * Type 2: Linear + 1D normal diffusion: Red.
%               * Type 3: Linear + 1D super diffusion: Green.
%               * Type 4: Linear + too short for 1D classification: Yellow.
%               * Type 5: Random/Unclassified + 2D confined diffusion: Blue.
%               * Type 6: Random/Unclassified + 2D normal diffusion: Cyan.
%               * Type 7: Random/Unclassified + 2D super diffusion: Magenta.
%               * Type 8: Random + too short for 2D classification: Purple.
%               * Type 0: Too short for any analysis: Light pink.
%
%Khuloud Jaqaman, August 2007

%% input - basic


%define colors to loop through
colorLoop = [1 0.7 0.7; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1]; %colors: 'light pink',r,g,b,y,m,c

%% store track positions, get track status and point status

schemeNumber=costMatLinkParam.linearMotion+1;

p = inputParser;
addOptional(p,'displayCutOff',1)
addOptional(p,'displayPred',1);
addOptional(p,'watchSchemeNumber',schemeNumber+1);        % Scheme to plot in the independant KF case.
addOptional(p,'showOneTrack',0);
parse(p,varargin{:});


% Here we put every needed data in a same shaped data structure in order to
% easily associate a feature with its scheme and coordinate
schemeAll = NaN(size(tracksFeatIndxLink));
posNoiseVarAll = NaN(size(tracksFeatIndxLink));
broNoiseVarAll = NaN(size(tracksFeatIndxLink));
xAll = NaN(size(tracksFeatIndxLink));
yAll = NaN(size(tracksFeatIndxLink));
dxAll = NaN(size(tracksFeatIndxLink));
dyAll = NaN(size(tracksFeatIndxLink));

startTime=startend(1);
endTime=startend(2);
probDim=2;
for t = startTime : endTime
    schemeAll((tracksFeatIndxLink(:,t)~=0),t) = KalmanInfoLink(t).scheme(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),2);
    if ndims(KalmanInfoLink(t).noiseVar)>3
        posNoiseVarAll((tracksFeatIndxLink(:,t)~=0),t) = ...
            KalmanInfoLink(t).noiseVar(1,1,tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),p.Results.watchSchemeNumber); 
        broNoiseVarAll((tracksFeatIndxLink(:,t)~=0),t) = ...
            KalmanInfoLink(t).noiseVar(1,1,tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),schemeNumber); 
        dxAll((tracksFeatIndxLink(:,t)~=0),t) = ...
            KalmanInfoLink(t).stateVec(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),3,1); 
        dyAll((tracksFeatIndxLink(:,t)~=0),t) = ...
            KalmanInfoLink(t).stateVec(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),4,1);         

    else
        posNoiseVarAll((tracksFeatIndxLink(:,t)~=0),t) = KalmanInfoLink(t).noiseVar(1,1,tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t)); 
        dxAll((tracksFeatIndxLink(:,t)~=0),t) = ...
            KalmanInfoLink(t).stateVec(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),3); 
        dyAll((tracksFeatIndxLink(:,t)~=0),t) = ...
            KalmanInfoLink(t).stateVec(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t),4);  
    end
    xAll((tracksFeatIndxLink(:,t)~=0),t) = movieInfo(t).xCoord(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t)); 
    yAll((tracksFeatIndxLink(:,t)~=0),t) = movieInfo(t).yCoord(tracksFeatIndxLink((tracksFeatIndxLink(:,t)~=0),t)); 

end

%copy brownStdMult into vector
kalmanStd = sqrt(probDim * abs(posNoiseVarAll));

stdMultInd = costMatLinkParam.brownStdMult*ones(size(posNoiseVarAll));

%get the search radius of each feature  and make sure it falls
%within reasonable limits
searchRadius = stdMultInd .* kalmanStd;
searchRadius((searchRadius>costMatLinkParam.maxSearchRadius)) = costMatLinkParam.maxSearchRadius;
searchRadius((searchRadius<costMatLinkParam.minSearchRadius)) = costMatLinkParam.minSearchRadius;

% % building segments data structure
% % data structure is build to fit nicely in plot wich 

%% make movie
numFramesMovie = diff(startend) + 1;
plotFullScreen=0;



if iscell(firstImageFile)
    [fpath,fname,fno,fext]=getFilenameBody(firstImageFile{1});
    dirName=[fpath,filesep];
    fName=[fname,fno,fext];
elseif ischar(firstImageFile)
    [fpath,fname,fno,fext]=getFilenameBody(firstImageFile);
    dirName=[fpath,filesep];
    fName=[fname,fno,fext];
end

movieType='mov';
%if input is valid ...
if(isa(fName,'char') && isa(dirName,'char'))
    
    %get all file names in stack
    outFileList = getFileStackNames([dirName,fName]);
    numFiles = length(outFileList);
    
    %determine     diffAnalysisRes = diffAnalysisRes(indx);which frames the files correspond to, and generate the inverse map
    %indicate missing frames with a zero
    frame2fileMap = zeros(numFiles,1);
    for iFile = 1 : numFiles
        [~,~,frameNumStr] = getFilenameBody(outFileList{iFile});
        frameNum = str2double(frameNumStr);
        frame2fileMap(frameNum) = iFile;
    end
    
    %assign as number of frames the last frame number observed
    numFrames = frameNum;
    
    %read first image to get image size
    currentImage = imread(outFileList{1});
    [isx,isy] = size(currentImage);
    
else %else, exit
    
    disp('--overlayTracksMovieNew: Bad file selection');
    return
    
end

[isx,isy] = size(currentImage);
if nargin < 10 || isempty(imageRange)
    imageRange = [1 isx; 1 isy];
end


xAll = xAll - (imageRange(2,1)-1);
yAll = yAll - (imageRange(1,1)-1);


intensityMinMax=[];

%initialize movie if it is to be saved
if saveMovie
    if movieType ~= 'png'
        movieVar = struct('cdata',[],'colormap',[]);
        movieVar = movieInfrastructure('initialize',movieType,dir2saveMovie,...
                                       movieName,numFramesMovie,movieVar,[]);
    end 
end



%keep only the frames of interest
outFileList = outFileList(frame2fileMap(startend(1)):frame2fileMap(startend(2)));
frame2fileMap = frame2fileMap(startend(1):startend(2));
indxNotZero = find(frame2fileMap~=0);
frame2fileMap(indxNotZero) = frame2fileMap(indxNotZero) - frame2fileMap(indxNotZero(1)) + 1;

%go over all specified frames
if plotFullScreen
    scrsz = get(0,'ScreenSize');
    h     = figure();
    set(h,'Position',scrsz);
else
    h = figure();
end


for iFrame = 1 : numFramesMovie
    
    if frame2fileMap(iFrame) ~= 0 %if frame exists
        
        %read specified image
        imageStack = imread(outFileList{frame2fileMap(iFrame)});
        
        %filter images if requested
        if filterSigma
            imageStack = filterGauss2D(imageStack,filterSigma);
        end
        
    else %otherwise
        
        %make empty frame
        imageStack = zeros(isx,isy);
        
    end
    
    %crop to region of interest
    imageStack = imageStack(imageRange(1,1):imageRange(1,2),...
        imageRange(2,1):imageRange(2,2),:);
    
    %     tmp = double(imageStack(:,:,1));
    %     minTmp = min(tmp(:));
    %     maxTmp = max(tmp(:));
    %     tmp = (tmp - minTmp)/(maxTmp - minTmp);
    %     imageStack(:,:,1) = uint8(tmp*255);
    %
    %     tmp = double(imageStack(:,:,2));
    %     minTmp = min(tmp(:));
    %     maxTmp = max(tmp(:));
    %     tmp = (tmp - minTmp)/(maxTmp - minTmp);
    %     imageStack(:,:,2) = uint8(tmp*255);
    %     imageStack(:,:,3) = imageStack(:,:,2);
    
    %plot image in current frame and show frame number
    clf;
    switch showRaw
        case 1
            
            axes('Position',[0 0 0.495 1]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            %             text(imageRange(2,1)+textDeltaCoord,imageRange(1,1)+...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            text(textDeltaCoord,...
                textDeltaCoord,num2str(iFrame+startend(1)-1),...
                'Color','white','FontSize',18);
            axes('Position',[0.505 0 0.495 1]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            
        case 2
            axes('Position',[0 0.505 1 0.495]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            %             text(imageRange(2,1)+textDeltaCoord,imageRange(1,1)+...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            text(textDeltaCoord,...
                textDeltaCoord,num2str(iFrame+startend(1)-1),...
                'Color','white','FontSize',18);
            %             text(textDeltaCoord-1,...
            %                 textDeltaCoord+2,[num2str((iFrame+startend(1)-2)*0.025,'%5.3f') ' s'],...
            %                 'Color','white','FontSize',18);
            axes('Position',[0 0 1 0.495]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
        otherwise
            axes('Position',[0 0 1 1]);
            imshow(imageStack,intensityMinMax);
            %             xlim(imageRange(2,:));
            %             ylim(imageRange(1,:));
            hold on;
            textDeltaCoord = min(diff(imageRange,[],2))/20;
            %             text(imageRange(2,1)+textDeltaCoord,imageRange(1,1)+...
            %                 textDeltaCoord,num2str(iFrame+startend(1)-1),...
            %                 'Color','white','FontSize',18);
            text(textDeltaCoord,...
                textDeltaCoord,num2str(iFrame+startend(1)-1),...
                'Color','white','FontSize',18);
    end
    
    %get tracks to plot
    plotOrNot = 0;
    if dragtailLength >= 0 %to plot tracks dynamically
        
        if iFrame > 1 
            % Extract segment to be plotted on current frame
            
            % building segment for plot input, keep only tracks
            % containing a spot on the current plotted frame.

            xStart=xAll(~isnan(xAll(:,iFrame)),max(1,iFrame-dragtailLength):(iFrame-1));
            yStart=yAll(~isnan(xAll(:,iFrame)),max(1,iFrame-dragtailLength):(iFrame-1));                        
            
            xEnd=xAll(~isnan(xAll(:,iFrame)),max(2,iFrame-dragtailLength+1):(iFrame));
            yEnd=yAll(~isnan(xAll(:,iFrame)),max(2,iFrame-dragtailLength+1):(iFrame));                  
            
            schemes=schemeAll(~isnan(xAll(:,iFrame)),max(1,iFrame-dragtailLength):(iFrame-1));           
                      
            %keep only feature with correspondance
            indxEndExist=~isnan(xEnd);
            indxStartExist=~isnan(xStart);  
            xStart = xStart(indxStartExist&indxEndExist)';
            xEnd = xEnd(indxStartExist&indxEndExist)';
            yStart = yStart(indxStartExist&indxEndExist)';
            yEnd = yEnd(indxStartExist&indxEndExist)';            
            schemes = schemes(indxStartExist&indxEndExist)';
            plotOrNot=1;
            
        end
        
    elseif dragtailLength == -1 %to plot tracks statically

    elseif dragtailLength == -2 %to plot tracks dynamically but keep them after they disappear

        
    end
    
    %plot tracks
    if plotOrNot
        
        %plot basic tracks
        plot([xStart(schemes==0);xEnd(schemes==0)],[yStart(schemes==0);yEnd(schemes==0)],...
         'Color',[1 0.7 0.7],'LineWidth',2); 
        if schemeNumber > 2
        plot([xStart(schemes==(schemeNumber-2));xEnd(schemes==(schemeNumber-2))],[yStart(schemes==(schemeNumber-2));yEnd(schemes==(schemeNumber-2))],...
         'Color',[1 0.7 0],'LineWidth',2); 
        plot([xStart(schemes==(2*schemeNumber-2));xEnd(schemes==(2*schemeNumber-2))],[yStart(schemes==(2*schemeNumber-2));yEnd(schemes==(2*schemeNumber-2))],...
             'Color','b','LineWidth',2);    
        end
        if schemeNumber > 1
            plot([xStart(schemes==(schemeNumber-1));xEnd(schemes==(schemeNumber-1))],[yStart(schemes==(schemeNumber-1));yEnd(schemes==(schemeNumber-1))],...
         'Color','r','LineWidth',2); 
            plot([xStart(schemes==(2*schemeNumber-1));xEnd(schemes==(2*schemeNumber-1))],[yStart(schemes==2*schemeNumber-1);yEnd(schemes==2*schemeNumber-1)],...
             'Color','y','LineWidth',2);    
        end 
        plot([xStart(schemes==schemeNumber);xEnd(schemes==schemeNumber)],[yStart(schemes==(schemeNumber));yEnd(schemes==schemeNumber)],...
             'Color','g','LineWidth',2);   
        plot([xStart(schemes==2*schemeNumber);xEnd(schemes==2*schemeNumber)],[yStart(schemes==2*schemeNumber);yEnd(schemes==2*schemeNumber)],...
             'Color',[ 1.0000    0.3333    0.6667],'LineWidth',2); % 
                        
        %light pink; the artificial repetition of the first line is for avoiding a mess in the first frame when tracks are not color-coded individually
%         plot([xStart;xEnd],yCoord2plot1,'Color',[1 0.7 0],'LineWidth',2); %orange
%         plot(xCoord2plot2,yCoord2plot2,'Color','r','LineWidth',2); %[1 0 0]
%         plot(xCoord2plot3,yCoord2plot3,'Color','g','LineWidth',2); %[0 1 0]

        
    end
        %plot points (features + gaps + merges + splits)
    if ~onlyTracks

        %red circles: detected feature in the middle of track with status 0
        if(p.Results.displayCutOff)
            if iFrame > 1
%                 if ndims(KalmanInfoLink(iFrame-1).noiseVar)==4
%                     noisevar=KalmanInfoLink(iFrame-1).noiseVar(1,1,:,p.Results.schemeNumber);
%                 else
%                     noisevar=KalmanInfoLink(iFrame-1).noiseVar(1,1,:);
%                 end 

                circles(searchRadius(:,iFrame-1), ...
                        [xAll(:,iFrame-1) yAll(:,iFrame-1)],'color','b');
                circles(3*sqrt(2*(broNoiseVarAll(:,iFrame-1))), ...
                        [xAll(:,iFrame-1) yAll(:,iFrame-1)],'color','g');
             end
        end 
        plot(xAll(:,iFrame),yAll(:,iFrame),'ro','MarkerSize',6);
        if (p.Results.displayPred)
            if iFrame>1
                if(p.Results.showOneTrack >0)
                    plot(xAll(p.Results.showOneTrack,iFrame-1)+dxAll(p.Results.showOneTrack,iFrame-1), ... 
                         yAll(p.Results.showOneTrack,iFrame-1)+dyAll(p.Results.showOneTrack,iFrame-1),'yo','MarkerSize',6);
                else 
                    plot(xAll(:,iFrame-1)+dxAll(:,iFrame-1),yAll(:,iFrame-1)+dyAll(:,iFrame-1),'yo','MarkerSize',6);
                end 
            end 
        end 

        %plot(xAll(points2plot,iFrame),yAll(points2plot,iFrame),'ro','MarkerSize',2);
    end
    
 
    
    %add frame to movie if movie is saved
    if saveMovie
        if movieType == 'tif'
            frameDir = [dir2saveMovie filesep movieName];
            [~,~] = mkdir(frameDir);
            ndigit = num2str(ceil(log10(numFramesMovie)));
            kformat = ['%.' ndigit 'd'];
            ext = '.tif';
            dest = [frameDir filesep 'frame_' num2str(iFrame, kformat) ext];
            hgexport(h, dest, hgexport('factorystyle'), 'Format', 'tiff');  
        else
        movieVar = movieInfrastructure('addFrame',movieType,dir2saveMovie,...
            movieName,numFramesMovie,movieVar,iFrame);
        end
    end
    
    %pause for a moment to see frame
    pause(0.01);
    
end %(for iFrame = 1 : numFramesMovie)

%finish movie
if saveMovie
    if movieType ~= 'png'
        movieInfrastructure('finalize',movieType,dir2saveMovie,...
                            movieName,numFramesMovie,movieVar,[]);
    end 
end

%% ~~~ end ~~~

