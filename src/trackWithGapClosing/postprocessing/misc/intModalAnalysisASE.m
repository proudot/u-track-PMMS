function [allIntensity,startIntensity,endIntensity] = intModalAnalysisASE(...
    tracks,alpha,variableMean,variableStd,maxNumGauss,plotRes,movieName)
%INTMODALANALYSISASE looks for intensity distribution modes at track starts and ends and throughout track lifetimes
%
%Synposis [allIntensity,startIntensity,endIntensity] = intModalAnalysisASE(...
%    tracks,alpha,variableMean,variableStd,maxNumGauss,plotRes,movieName)
%
%INPUT  tracks      : Output of trackCloseGapsKalman
%       alpha       : Alpha-value for the statistical test that compares the
%                     fit of n+1 Gaussians to the fit of n Gaussians. Can
%                     either contain two values, one for fitting the
%                     intensity distribution, one for fitting the
%                     distribution of Gaussian means, or one value, which
%                     will then be used for both.
%                     Optional. Default: [0.05 0.05].
%       variableMean: 0 if assuming the fixed relationship
%                     (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                     1 if there is no relationship between the means of
%                     different Gaussians. Can either contain two values,
%                     one for fitting the intensity distribution, one for
%                     fitting the distribution of Gaussian means, or one
%                     value, which will then be used for both.
%                     Optional. Default: [0 0].
%       variableStd : 0 if assuming that all Gaussians have the same
%                     standard deviation. 1 if there is no relationship 
%                     between the standard deviations of different
%                     Gaussians, 2 if assuming the relationship 
%                     (std of nth Gaussian) = sqrt(n) * (std of 1st Gaussian). 
%                     variableStd can equal 2 only if variableMean is 0.
%                     Can either contain two values, one for fitting the
%                     intensity distribution, one for fitting the
%                     distribution of Gaussian means, or one value, which
%                     will then be used for both.
%                     Optional. Default: [0 0].
%       maxNumGauss : Upper limit to the number of fitted Gaussians.
%                     Optional. Default: 10.
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0.
%       movieName   : Name of movie. Needed to put as title of plot if plotting.
%                     Optional. Default: 'movie1'.
%
%OUTPUT allIntensity, startIntensity, endIntensity: Structures with
%                 analysis results for intensities throughout the whole
%                 movie, just after tracks start, and just before tracks
%                 end, respectively, containing the fields:
%           .numGauss: Number of Gaussians best fitting the intensity
%                      distribution.
%           .firstGaussMean: Mean of the first Gaussian fitting the
%                      intensity distribution.
%           .firstGaussStd: Standard deviation of the first Gaussian.
%           .ratioPeak2toPeak1: Ratio of amplitude of second Gaussian to
%                      first Gaussian. NaN indicates that a frame does not
%                      have a second Gaussian.
%           .gaussFit.param: Parameters of the Gaussian fit for each intensity
%                      histogram. Equivalent to the output "gaussParam" of
%                      fitHistWithGaussians.
%           .gaussMeanModes: Parameters of the Gaussian fit to the
%                      distribution of Gaussian means. Equivalent to the
%                      output "gaussParam" of fitHistWithGaussians.
%
%Khuloud Jaqaman, September 2009

%% Output

allIntensity = [];
startIntensity = [];
endIntensity = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--intModalAnalysisASE: Function needs at least 1 input argument!');
    return
end

%check alpha
if nargin < 2 || isempty(alpha)
    alpha = [0.05 0.05];
elseif length(alpha) == 1
    alpha = alpha * ones(1,2);
end

%check variableMean
if nargin < 3 || isempty(variableMean)
    variableMean = [0 0];
elseif length(variableMean) == 1
    variableMean = variableMean * ones(1,2);
end

%check variableStd
if nargin < 4 || isempty(variableStd)
    variableStd = [0 0];
elseif length(variableStd) == 1
    variableStd = variableStd * ones(1,2);
end

%check maxNumGauss
if nargin < 5 || isempty(maxNumGauss)
    maxNumGauss = 10;
end

%check plotRes
if nargin < 6 || isempty(plotRes)
    plotRes = 0;
end

%check movieName
if nargin < 7 || isempty(movieName)
    movieName = 'movie 1';
end    

%define number of subsets for histogram fitting
numSubsets = 500;

%definesubset size for histogram fitting
subsetSize = 200;

%% Preamble

%find the length of the movie
seqOfEventsAll = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEventsAll(:,1));

%find tracks with only 1 segment and that start and end within the movie
criteria.numSegments.max = 1;
criteria.startTime.min = 2;
criteria.endTime.max = numFrames - 1;
indxOneSegment = chooseTracks(tracks,criteria);
numTracksOneSegment = length(indxOneSegment);

%convert these tracks into matrix format and extract intensities
tracksOneSegmentMat = convStruct2MatNoMS(tracks(indxOneSegment));
ampOneSegmentMat = tracksOneSegmentMat(:,4:8:end);

%put all the intensities in a vector with no NaNs
ampAllVec = ampOneSegmentMat(:);
ampAllVec = ampAllVec(~isnan(ampAllVec));
numAmpAll = length(ampAllVec);

%find the start times and end times of these tracks
trackSEL = getTrackSEL(tracksOneSegmentMat);
trackStartTime = trackSEL(:,1);
trackEndTime = trackSEL(:,2);

%get the track intensities just after the tracks start and just before they
%end
ampStart = NaN(numTracksOneSegment,1);
ampEnd = NaN(numTracksOneSegment,1);
for iTrack = 1 : numTracksOneSegment
    ampStart(iTrack) = ampOneSegmentMat(iTrack,trackStartTime(iTrack));
    ampEnd(iTrack) = ampOneSegmentMat(iTrack,trackEndTime(iTrack));
end

% %convert all tracks into matrix format and extract intensities
% tracksAllMat = convStruct2MatIgnoreMS(tracks);
% ampAllVec = tracksAllMat(:,4:8:end);
% ampAllVec = ampAllVec(:);
% ampAllVec = ampAllVec(~isnan(ampAllVec));
% numAmpAll = length(ampAllVec);

%% Overall distribution analysis

%if sample size is smaller than subset size, don't choose subsets
%just fit once and that's all
if numAmpAll <= subsetSize
    numSubsetsTmpA = 1;
    subsetSizeTmpA = numAmpAll;
else
    numSubsetsTmpA = numSubsets;
    subsetSizeTmpA = subsetSize;
end
    
%choose subsets of the amplitudes and try to fit multiple
%Gaussians to the intensity histogram
[numGauss,firstGaussMean,firstGaussStd,ratioPeak2toPeak1] = deal(NaN(numSubsetsTmpA,1));
gaussFit = repmat(struct('param',[]),numSubsetsTmpA,1);
for iSubset = 1 : numSubsetsTmpA
    indx = randsample(numAmpAll,subsetSizeTmpA);
    [dummy,dummy,gaussParamSubset] = fitHistWithGaussians(...
        ampAllVec(indx),alpha(1),variableMean(1),variableStd(1),0,maxNumGauss,2);
    numGauss(iSubset) = size(gaussParamSubset,1);
    firstGaussMean(iSubset) = gaussParamSubset(1,1);
    firstGaussStd(iSubset) = gaussParamSubset(1,2);
    if numGauss(iSubset) > 1
        ratioPeak2toPeak1(iSubset) = ...
            (gaussParamSubset(2,3)/gaussParamSubset(2,2))/...
            (gaussParamSubset(1,3)/gaussParamSubset(1,2));
    end
    gaussFit(iSubset).param = gaussParamSubset;
end

%collect all the Gaussian means into a distribution
gaussMeanDistr = catStruct(1,'gaussFit.param(:,1)');

%perform modal analysis on the distribution of the means
[dummy,dummy,gaussParamSubset] = fitHistWithGaussians(gaussMeanDistr,...
    alpha(2),variableMean(2),variableStd(2),plotRes,maxNumGauss,2,[movieName ...
    ' - Modal analysis of Gaussian means, all intensities']);

%prepare output
allIntensity.numGauss = numGauss;
allIntensity.firstGaussMean = firstGaussMean;
allIntensity.firstGaussStd = firstGaussStd;
allIntensity.ratioPeak2toPeak1 = ratioPeak2toPeak1;
allIntensity.gaussFit = gaussFit;
allIntensity.gaussMeanModes = gaussParamSubset;

%% Start distribution analysis

%if sample size is smaller than subset size, don't choose subsets
%just fit once and that's all
if numTracksOneSegment <= subsetSize
    numSubsetsTmpS = 1;
    subsetSizeTmpS = numTracksOneSegment;
else
    numSubsetsTmpS = numSubsets;
    subsetSizeTmpS = subsetSize;
end

%choose subsets of these amplitudes and try to fit multiple
%Gaussians to the intensity histogram
[numGauss,firstGaussMean,firstGaussStd,ratioPeak2toPeak1] = deal(NaN(numSubsetsTmpS,1));
gaussFit = repmat(struct('param',[]),numSubsetsTmpS,1);
for iSubset = 1 : numSubsetsTmpS
    indx = randsample(numTracksOneSegment,subsetSizeTmpS);
    [dummy,dummy,gaussParamSubset] = fitHistWithGaussians(...
        ampStart(indx),alpha(1),variableMean(1),variableStd(1),0,maxNumGauss,2);
    numGauss(iSubset) = size(gaussParamSubset,1);
    firstGaussMean(iSubset) = gaussParamSubset(1,1);
    firstGaussStd(iSubset) = gaussParamSubset(1,2);
    if numGauss(iSubset) > 1
        ratioPeak2toPeak1(iSubset) = ...
            (gaussParamSubset(2,3)/gaussParamSubset(2,2))/...
            (gaussParamSubset(1,3)/gaussParamSubset(1,2));
    end
    gaussFit(iSubset).param = gaussParamSubset;
end

%collect all the Gaussian means into a distribution
gaussMeanDistr = catStruct(1,'gaussFit.param(:,1)');

%perform modal analysis on the distribution of the means
[dummy,dummy,gaussParamSubset] = fitHistWithGaussians(gaussMeanDistr,...
    alpha(2),variableMean(2),variableStd(2),plotRes,maxNumGauss,2,[movieName ...
    ' - Modal analysis of Gaussian means, start intensities']);

%prepare output
startIntensity.numGauss = numGauss;
startIntensity.firstGaussMean = firstGaussMean;
startIntensity.firstGaussStd = firstGaussStd;
startIntensity.ratioPeak2toPeak1 = ratioPeak2toPeak1;
startIntensity.gaussFit = gaussFit;
startIntensity.gaussMeanModes = gaussParamSubset;

%% End distribution analysis

%if sample size is smaller than subset size, don't choose subsets
%just fit once and that's all
if numTracksOneSegment <= subsetSize
    numSubsetsTmpE = 1;
    subsetSizeTmpE = numTracksOneSegment;
else
    numSubsetsTmpE = numSubsets;
    subsetSizeTmpE = subsetSize;
end

%choose subsets of these amplitudes and try to fit multiple
%Gaussians to the intensity histogram
[numGauss,firstGaussMean,firstGaussStd,ratioPeak2toPeak1] = deal(NaN(numSubsetsTmpE,1));
gaussFit = repmat(struct('param',[]),numSubsetsTmpE,1);
for iSubset = 1 : numSubsetsTmpE
    indx = randsample(numTracksOneSegment,subsetSizeTmpE);
    [dummy,dummy,gaussParamSubset] = fitHistWithGaussians(...
        ampEnd(indx),alpha(1),variableMean(1),variableStd(1),0,maxNumGauss,2);
    numGauss(iSubset) = size(gaussParamSubset,1);
    firstGaussMean(iSubset) = gaussParamSubset(1,1);
    firstGaussStd(iSubset) = gaussParamSubset(1,2);
    if numGauss(iSubset) > 1
        ratioPeak2toPeak1(iSubset) = ...
            (gaussParamSubset(2,3)/gaussParamSubset(2,2))/...
            (gaussParamSubset(1,3)/gaussParamSubset(1,2));
    end
    gaussFit(iSubset).param = gaussParamSubset;
end

%collect all the Gaussian means into a distribution
gaussMeanDistr = catStruct(1,'gaussFit.param(:,1)');

%perform modal analysis on the distribution of the means
[dummy,dummy,gaussParamSubset] = fitHistWithGaussians(gaussMeanDistr,...
    alpha(2),variableMean(2),variableStd(2),plotRes,maxNumGauss,2,[movieName ...
    ' - Modal analysis of Gaussian means, end intensities']);

%prepare output
endIntensity.numGauss = numGauss;
endIntensity.firstGaussMean = firstGaussMean;
endIntensity.firstGaussStd = firstGaussStd;
endIntensity.ratioPeak2toPeak1 = ratioPeak2toPeak1;
endIntensity.gaussFit = gaussFit;
endIntensity.gaussMeanModes = gaussParamSubset;

%% Plotting

if plotRes
    
    %start new figure and put title
    figure('Name',movieName,'NumberTitle','off');
    
    %get maximum number of fitted Gaussians
    maxNumGauss = max([allIntensity.numGauss; startIntensity.numGauss; ...
        endIntensity.numGauss]);
    
    %get maximum Gaussian mean
    gaussMeanDistr = [catStruct(1,'allIntensity.gaussFit.param(:,1)'); ...
        catStruct(1,'startIntensity.gaussFit.param(:,1)'); ...
        catStruct(1,'endIntensity.gaussFit.param(:,1)')];
    maxGaussMean = max(gaussMeanDistr);
    
    %plot number of fitted Gaussians for all intensities
    subplot(2,3,1)
    hold on
    plot(allIntensity.numGauss)
    axis([1 numSubsetsTmpA 0 maxNumGauss+1])
    xlabel('Subset number')
    ylabel('Number of Modes')
    title('All intensities')
    
    %plot number of fitted Gaussians for start intensities
    subplot(2,3,2)
    hold on
    plot(startIntensity.numGauss)
    axis([1 numSubsetsTmpS 0 maxNumGauss+1])
    xlabel('Subset number')
    ylabel('Number of Modes')
    title('Start intensities')

    %plot number of fitted Gaussians for end intensities
    subplot(2,3,3)
    hold on
    plot(endIntensity.numGauss)
    axis([1 numSubsetsTmpE 0 maxNumGauss+1])
    xlabel('Subset number')
    ylabel('Number of Modes')
    title('End intensities')
    
    %plot the Gaussian means for all intensities
    subplot(2,3,4)
    hold on
    for i=1:length(allIntensity.numGauss)
        plot(i*ones(allIntensity.numGauss(i),1),allIntensity.gaussFit(i).param(:,1),'.')
    end
    
    %plot on top of that the means of the Gaussian modes
    numModes = size(allIntensity.gaussMeanModes,1);
    xValues = repmat([1 length(allIntensity.numGauss)]',1,numModes);
    yValues = repmat(allIntensity.gaussMeanModes(:,1)',2,1);
    plot(xValues,yValues,'r')
    
    %add titles
    axis([1 numSubsetsTmpA 0 1.1*maxGaussMean])
    xlabel('Subset number')
    ylabel('Mode center')

    %plot the Gaussian means for start intensities
    subplot(2,3,5)
    hold on
    for i=1:length(startIntensity.numGauss)
        plot(i*ones(startIntensity.numGauss(i),1),startIntensity.gaussFit(i).param(:,1),'.')
    end
    
    %plot on top of that the means of the Gaussian modes
    numModes = size(startIntensity.gaussMeanModes,1);
    xValues = repmat([1 length(startIntensity.numGauss)]',1,numModes);
    yValues = repmat(startIntensity.gaussMeanModes(:,1)',2,1);
    plot(xValues,yValues,'r')
    
    %add titles
    axis([1 numSubsetsTmpS 0 1.1*maxGaussMean])
    xlabel('Subset number')
    ylabel('Mode center')

    %plot the Gaussian means for end intensities
    subplot(2,3,6)
    hold on
    for i=1:length(endIntensity.numGauss)
        plot(i*ones(endIntensity.numGauss(i),1),endIntensity.gaussFit(i).param(:,1),'.')
    end
    
    %plot on top of that the means of the Gaussian modes
    numModes = size(endIntensity.gaussMeanModes,1);
    xValues = repmat([1 length(endIntensity.numGauss)]',1,numModes);
    yValues = repmat(endIntensity.gaussMeanModes(:,1)',2,1);
    plot(xValues,yValues,'r')
    
    %add titles
    axis([1 numSubsetsTmpE 0 1.1*maxGaussMean])
    xlabel('Subset number')
    ylabel('Mode center')

end


% if plotRes
%     
%     %start new figure and put title
%     figure('Name',movieName,'NumberTitle','off');
%     
%     %plot distribution of number of fitted Gaussians
%     subplot(1,3,1)
%     hold on;
%     maxNumGauss = max([allIntensity.numGauss; endIntensity.numGauss]);
%     nAll = hist(allIntensity.numGauss,1:maxNumGauss);
%     nAll = nAll / sum(nAll);
%     nEnd = hist(endIntensity.numGauss,1:maxNumGauss);
%     nEnd = nEnd / sum(nEnd);
%     bar([nAll; nEnd]')
%     legend('All','End')
%     xlabel('Number of Gaussians')
%     ylabel('Normalized frequency')
%     
%     %plot the means and stds of the first Gaussian
%     subplot(1,3,2)
%     hold on
%     if numSubsetsTmp2 > 1
%         plot(allIntensity.firstGaussMean,'k')
%         plot(allIntensity.firstGaussStd,'r')
%     else
%         plot(allIntensity.firstGaussMean,'k+')
%         plot(allIntensity.firstGaussStd,'r+')
%     end
%     if numSubsetsTmp1 > 1
%         plot(endIntensity.firstGaussMean,'c')
%         plot(endIntensity.firstGaussStd,'m')
%     else
%         plot(endIntensity.firstGaussMean,'c+')
%         plot(endIntensity.firstGaussStd,'m+')
%     end
%     legend('Mean, all','Std, all','Mean, end','Std, end')
%     numSubsetsMax = max([numSubsetsTmp1 numSubsetsTmp2]);
%     plot([0 numSubsetsMax],mean(allIntensity.firstGaussMean)*[1 1],'k')
%     plot([0 numSubsetsMax],mean(allIntensity.firstGaussStd)*[1 1],'r')
%     plot([0 numSubsetsMax],mean(endIntensity.firstGaussMean)*[1 1],'c')
%     plot([0 numSubsetsMax],mean(endIntensity.firstGaussStd)*[1 1],'m')
%     xlabel('Subset number')
%     ylabel('Mean/std of first Gaussian')
%     
%     %plot the raio of the 2-fluorophore population to the 1-fluorophore
%     %population
%     subplot(1,3,3)
%     hold on
%     if numSubsetsTmp2 > 1
%         plot(allIntensity.ratioPeak2toPeak1,'g')
%     else
%         plot(allIntensity.ratioPeak2toPeak1,'g+')
%     end
%     if numSubsetsTmp1 > 1
%         plot(endIntensity.ratioPeak2toPeak1,'b')
%     else
%         plot(endIntensity.ratioPeak2toPeak1,'b+')
%     end
%     legend('All','End')
%     plot([0 numSubsetsMax],nanmean(allIntensity.ratioPeak2toPeak1)*[1 1],'g')
%     plot([0 numSubsetsMax],nanmean(endIntensity.ratioPeak2toPeak1)*[1 1],'b')
%     xlabel('Subset number')
%     ylabel('2-fluor. : 1-fluor. ratio')
%     
% end


