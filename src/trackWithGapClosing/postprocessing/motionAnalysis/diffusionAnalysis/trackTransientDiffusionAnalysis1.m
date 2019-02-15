function [transDiffAnalysisRes,errFlag] = trackTransientDiffusionAnalysis1(tracks,...
    extractType,probDim,alphaValues,minDuration,plotRes,confRadMin)
%TRACKTRANSIENTDIFFUSIONANALYSIS1 performs a rolling window MSS diffusion analysis
%
%SYNOPSIS [transDiffAnalysisRes,errFlag] = trackTransientDiffusionAnalysis1(tracks,...
%     extractType,probDim,alphaValues,minDuration,plotRes,confRadMin)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits).
%       extractType : 1 - Analyze every track segment separately.
%                     2 - Extract from each compound track the longest
%                         trajectory to use in analysis - NOT IMPLEMENTED
%                         YET.
%                     Variable irrelevant if tracks are input as a matrix.
%                     Optional. Default: 1.
%       probDim     : Problem dimensionality.
%                     Optional. Default: 2.
%       alphaValues : Alpha-value for classification. Can take the values
%                     0.2, 0.1, 0.05 and 0.01. One can enter one value, in
%                     which case it will be used for both confined and
%                     directed, or two values, where the first will be used
%                     for confined and the second for directed.
%                     Optional. Default: 0.05 for both.
%       minDuration : Minimum category duration. One can enter one value,
%                     in which case it will be used for both confined and
%                     directed, or two values, where the first will be used
%                     for confined and the second for directed.
%                     Optional. Default: [8 2].
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0.
%                     Results can be plotted only if problem is 2D.
%                     color-coding:
%                     *blue: confined diffusion.
%                     *cyan: normal diffusion.
%                     *magenta: super diffusion.
%                     *black: unclassified.
%       confRadMin  : 1 to calculate the confinement radius of confined
%                     particles using the minimum positional standard
%                     deviation, 0 to calculate it using the mean
%                     positional standard deviation.
%                     Optional. Default: 0.
%
%OUTPUT transDiffAnalysisRes : And array (size = number of tracks) with the
%                     field ".segmentClass", which contains the fields:
%           .momentScalingSpectrum: (Number of classification
%                     subparts)-by-(20+probDim) matrix, where each row
%                     contains the information for one classification
%                     subpart, and the columns store the following:
%                     (1) Start frame of subpart.
%                     (2) End frame of subpart.
%                     (3) Classification of subpart: 1 = confined, 2 =
%                         free, 3 = directed.
%                     (4) MSS slope resulting in classification.
%                     (5-11) Generalized diffusion coefficients for moment
%                            orders 0-6.
%                     (12-18) Scaling power for moment orders 0-6.
%                     (19) Normal diffusion coefficient (from the MSD).
%                     (20) Confinement radius, if subpart is classified as
%                          confined (NaN otherwise).
%                     (21/22/23) Center of subpart, if subpart is
%                                classified as confined (NaN otherwise).
%           .asymmetry, .asymmetryAfterMSS: NOT IMPLEMENTED RIGHT NOW.
%
%       errFlag         : 0 if executed normally, 1 otherwise.
%
%REMARKS While tracks do not have to be linear in order to be asymmetric,
%the last analysis step assumes that tracks are linear.
%
%Khuloud Jaqaman, March 2008

%% Output

transDiffAnalysisRes = [];
errFlag = 0;

%% Input

%check whether tracks were input
if nargin < 1
    disp('--trackTransientDiffusionAnalysis1: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--trackTransientDiffusionAnalysis1: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(alphaValues)
    alphaValues = [0.05 0.05];
elseif length(alphaValues) == 1
    alphaValues = [alphaValues alphaValues];
end

if nargin < 5 || isempty(minDuration)
    minDuration = [8 2];
elseif length(minDuration) == 1
    minDuration = [minDuration minDuration];
end

if nargin < 6 || isempty(plotRes)
    plotRes = 0;
elseif plotRes == 1 && probDim ~= 2
    disp('--trackTransientDiffusionAnalysis1: Cannot plot tracks if problem is not 2D!');
    plotRes = 0;
end

if nargin < 7 || isempty(confRadMin)
    confRadMin = 0;
end

if errFlag
    disp('--trackTransientDiffusionAnalysis1: Please fix input variables');
    return
end

%THE ASYMMETRY ANALYSIS PART NEEDS WORK, SO I WILL IMPOSE NOT USING IT FOR
%NOW.
checkAsym = 0;

%define window sizes
windowAsym = 5;
windowMSS = 21;
windowMSSMin = 20;
% halfWindowAsym = (windowAsym - 1) / 2;
halfWindowMSS = (windowMSS - 1) / 2;

%define duration of a Brownian segment that leads to lumping it with the
%segments around it
windowBrown = 20;

%specify MSS analysis moment orders
momentOrders = 0 : 6;

%% Track extraction for analysis

%store input tracks in a new variable
tracksInput = tracks;

%extract segments for analysis if tracks were input as a structure that
%might contain merges and splits
%the point is to reduce compound tracks that contain merges and splits into
%simple separate tracks
%thus this step is not necessary if the tracks were input as a matrix,
%which by definition does not contain unresolved compound tracks.
if isstruct(tracks)

    %get number of input tracks from structure
    numInputTracks = length(tracksInput);

    clear tracks

    switch extractType

        case 1 %retrieve every track segment separately

            [tracks,dummy,compTrackStartRow,numSegments] = ...
                convStruct2MatIgnoreMS(tracksInput);

        case 2 %make the longest track possible, given all the merges and splits

            disp('Sorry - not implemented yet!')
            errFlag = 1;
            return

    end

else

    %get number of input tracks from matrix
    numInputTracks = size(tracksInput,1);

    %indicate rows where tracks start (trivial in this case)
    compTrackStartRow = (1 : numInputTracks)';

    %indicate number of segments in each track (1 for all tracks)
    numSegments = ones(numInputTracks,1);

end

%get number of track segments to be analyzed
numTrackSegments = size(tracks,1);

%get track segment start times, end times and life times
trackSEL = getTrackSEL(tracks);

%find track segments that are long enough for analysis
if checkAsym
    indx4analysis = find(trackSEL(:,3) >= windowAsym);
else
    indx4analysis = find(trackSEL(:,3) >= windowMSSMin);
end
indxNot4analysis = setdiff((1:numTrackSegments)',indx4analysis);

%% Rolling window classification

%reserve memory
trackSegmentClassRes = repmat(struct('asymmetry',NaN(1,3),...
    'momentScalingSpectrum',NaN(1,21+probDim),...
    'asymmetryAfterMSS',NaN(1,3)),...
    numTrackSegments,1);

%go over all analyzable track segments
for iTrack = indx4analysis'

    %% Asymmetry analysis to get directed parts of the track segment

    %get track segment start, end and life times
    trackSELCurrent = trackSEL(iTrack,:);

    if checkAsym

        %         %get number of asymmetry analysis rolling windows for this track
        %         %segment
        %         numRollWindows = trackSELCurrent(3) - windowAsym + 1;
        %
        %         %initialize point classification vector
        %         pointClassAsym = NaN(trackSELCurrent(3),1);
        %
        %         %get the particle positions along the track segment
        %         coordX = tracks(iTrack,1:8:end)';
        %         coordY = tracks(iTrack,2:8:end)';
        %         coordZ = tracks(iTrack,3:8:end)';
        %         coordXYZ = [coordX coordY coordZ];
        %         coordXYZ = coordXYZ(trackSELCurrent(1):trackSELCurrent(2),:);
        %
        %         %go over all windows
        %         for iWindow = 1 : numRollWindows
        %
        %             %determine whether the track segment is sufficiently asymmetric
        %             %in this window
        %             coord2classify = coordXYZ(iWindow:iWindow+windowAsym-1,1:probDim);
        %             if any(~isnan(coord2classify))
        %                 [asymParamT,asymFlag] = asymDeterm2D3D(coord2classify,alphaAsym);
        %             else
        %                 asymFlag = NaN;
        %             end
        %
        %             %classify window middle point as ...
        %             %1 = linear, if the asymmetry parameter is larger than the threshold
        %             %0 = not linear, if the asymmetry parameter is smaller than the
        %             %threshold
        %             %otherwise, keep track segment classification as undetermined
        %             pointClassAsym(iWindow+(windowAsym-1)/2) = asymFlag;
        %
        %         end
        %
        %         %for any point classified as asymmetric, expand its classification
        %         %to the points surrounding it and corresponding to its window
        %         indxGood = find(pointClassAsym==1);
        %         for iPoint = indxGood'
        %             pointClassAsym(iPoint-(windowAsym-1)/2:iPoint+(windowAsym-1)/2) = 1;
        %         end
        %
        %         %expand the classification to the first and last few points of a
        %         %track segment
        %         pointClassAsym(1:(windowAsym-1)/2) = pointClassAsym((windowAsym+1)/2);
        %         pointClassAsym(end-(windowAsym-3)/2:end) = pointClassAsym(end-(windowAsym-1)/2);
        %
        %         %give unclassified points the classification of the point before them
        %         unclassPoints = find(isnan(pointClassAsym));
        %         if ~isempty(unclassPoints)
        %             if unclassPoints(1) == 1
        %                 unclassPoints = unclassPoints(2:end);
        %             end
        %             for iPoint = unclassPoints'
        %                 pointClassAsym(iPoint) = pointClassAsym(iPoint-1);
        %             end
        %         end
        %
        %         %if this track segment got some classification ...
        %         indxGood = find(~isnan(pointClassAsym));
        %         if ~isempty(indxGood)
        %
        %             %if there are unclassified points in the beginning, give them the
        %             %classification of the first classified point
        %             if isnan(pointClassAsym(1))
        %                 pointClassAsym(1:indxGood(1)-1) = pointClassAsym(indxGood(1));
        %             end
        %
        %             %take the pointClassAsym difference: 0 means staying in the same
        %             %state, 1 means switching from symmetric to asymmetric, and -1
        %             %means switching from asymmetric to symmetric
        %             %add one point to the beginning and the end, to get end points
        %             pointClassAsymDiff = diff([1-pointClassAsym(1); pointClassAsym; 1-pointClassAsym(end)]);
        %             switch0110 = find(pointClassAsymDiff ~= 0);
        %
        %             %make an array of classes and their starting points
        %             %note that the very last class is not real
        %             switch0110 = [switch0110 [pointClassAsym(switch0110(1:end-1)); 1-pointClassAsym(end)]]; %#ok
        %
        %             %get the duration of each class
        %             classDuration = [diff(switch0110(:,1)) switch0110(1:end-1,2)];
        %
        %             %find all track segment parts originally classified as asymmetric but
        %             %for less than minDuration(2) frames, and convert them into symmetric
        %             badClass = find(classDuration(:,1) < minDuration(2) & classDuration(:,2)==1);
        %             for iClass = badClass'
        %                 pointClassAsym(switch0110(iClass,1):switch0110(iClass+1,1)-1) = 0;
        %             end
        %
        %             %take the pointClassAsym difference again and get class durations
        %             %again ...
        %             pointClassAsymDiff = diff([1-pointClassAsym(1); pointClassAsym; 1-pointClassAsym(end)]);
        %             switch0110 = find(pointClassAsymDiff ~= 0);
        %             switch0110 = [switch0110 [pointClassAsym(switch0110(1:end-1)); 1-pointClassAsym(end)]]; %#ok
        %             classDuration = [diff(switch0110(:,1)) switch0110(1:end-1,2)];
        %
        %             %find all track segment parts originally classified as
        %             %symmetric, but for less than windowAsym frames, and convert
        %             %them into asymmetric
        %             badClass = find(classDuration(:,1) < windowAsym & classDuration(:,2)==0);
        %             for iClass = badClass'
        %                 pointClassAsym(switch0110(iClass,1):switch0110(iClass+1)-1,1) = 1;
        %             end
        %
        %             %for the last time, find switching points
        %             %this is the final asymmetry classification of this track
        %             %segment
        %             pointClassAsymDiff = diff([1-pointClassAsym(1); pointClassAsym; 1-pointClassAsym(end)]);
        %             switch0110 = find(pointClassAsymDiff ~= 0);
        %             trackClassAsym = [switch0110(1:end-1) switch0110(2:end)-1 pointClassAsym(switch0110(1:end-1))];
        %
        %             %shift time to start track segment classification at track
        %             %segment start time
        %             trackClassAsym(:,1:2) = trackClassAsym(:,1:2) + trackSELCurrent(1) - 1;
        %
        %         else
        %
        %             %if there are no classified points, then keep whole track
        %             %segment as unclassified
        %             trackClassAsym = [trackSELCurrent(1:2) NaN];
        %
        %         end %(if ~isempty(indxGood) ... else ...)

    else

        %if no asymmetry classification is to be performed, keep whole
        %track segment as unclassified
        trackClassAsym = [trackSELCurrent(1:2) NaN];

    end %(if checkAsym ... else ...)

    %% MSS analysis to divide track segment into parts with different diffusion behavior

    %assign initial MSS classification
    %all parts classified above as asymmetric get a -1
    trackClassMSS = trackClassAsym;
    trackClassMSS(:,3) = -trackClassMSS(:,3);
    trackClassMSS(trackClassMSS(:,3)~=-1,3) = NaN;
    trackClassMSS(:,4:20+probDim) = NaN;
    oldPart2newPartMap = (1 : size(trackClassMSS,1))';

    %initialize expanded asymmetry classification matrix
    trackClassAsymExp = trackClassAsym;

    %find the length of each part
    trackPartLength = trackClassMSS(:,2) - trackClassMSS(:,1) + 1;

    %find track segment parts that are not classified as asymmetric and
    %that are longer than windowMSSMin
    trackParts2analyze = find(trackClassMSS(:,3) ~= -1 & trackPartLength >= windowMSSMin);

    %if there are parts to analyze ...
    if ~isempty(trackParts2analyze)

        %go over these parts ...
        for iPart = trackParts2analyze'

            %get starting point of this part
            trackPartStart = trackClassMSS(oldPart2newPartMap(iPart),1);

            %get number of MSS analysis rolling windows in this part
            numRollWindows = trackPartLength(iPart) - windowMSS + 1;

            %if number of rolling windows is larger than the minimum
            %required duration of a classification, proceed with rolling
            %window analysis
            if numRollWindows > minDuration(1)

                %initialize point classification vector
                pointClassMSS = NaN(numRollWindows,1);
                mssSlope = pointClassMSS;

                %go over all windows
                for iWindow = 1 : numRollWindows

                    %get window start and end points
                    startPoint = trackPartStart + iWindow - 1;
                    endPoint = startPoint + windowMSS - 1;

                    %call MSS analysis code
                    %classify window middle point as ...
                    %1 = confined diffusion, if MSS slope is smaller than lower threshold
                    %3 = directed, if MSS slope is larger than upper threshold
                    %2 = pure Brownian, if MSS slope is in between the two thresholds
                    %note that the first and last halfWindowMSS points are
                    %not included in pointClassMSS and mssSlope
                    [pointClassMSS(iWindow),mssSlope(iWindow)] = trackMSSAnalysis(...
                        tracks(iTrack,8*(startPoint-1)+1:8*endPoint),...
                        probDim,momentOrders,alphaValues);

                end %(for iWindow = 1 : numRollWindows)

                %replace NaN with 0 in pointClassMSS
                pointClassMSS(isnan(pointClassMSS)) = 0;

                %if this part got some classification ...
                if ~isempty(pointClassMSS~=0)
                    
                    %reclassify unclassified subparts
                    pointClassMSS = reclassUnclassPoints(pointClassMSS,...
                        mssSlope,tracks(iTrack,:),halfWindowMSS,probDim,...
                        momentOrders,alphaValues,1,trackPartStart);

                    %make sure that subparts classified as confined or
                    %directed span the minimum length
                    %if they do, expand them
                    [pointClassMSS,mssSlope,halfWindowStartRemoved,...
                        halfWindowEndRemoved] = reclassConfinedDirectedPoints(...
                        pointClassMSS,mssSlope,minDuration,halfWindowMSS);

                    %make sure that subparts classified as Brownian span
                    %the minimum length
                    pointClassMSS = reclassBrownianPoints(pointClassMSS,...
                        mssSlope,windowBrown);

                    %now add back the first and last halfWindowMSS points
                    %and give them the classification of the immediate
                    %neighboring points, if they haven't been classified
                    %yet
                    if halfWindowStartRemoved
                        pointClassMSS = [pointClassMSS(1)*ones(halfWindowMSS,1); ...
                            pointClassMSS]; %#ok<AGROW>
                    end
                    if halfWindowEndRemoved
                        pointClassMSS = [pointClassMSS; pointClassMSS(end)*...
                            ones(halfWindowMSS,1)]; %#ok<AGROW>
                    end

                    %try to reclassify any remaining unclassified subparts
                    %for the last time
                    pointClassMSS = reclassUnclassPoints(pointClassMSS,...
                        mssSlope,tracks(iTrack,:),halfWindowMSS,probDim,...
                        momentOrders,alphaValues,0,trackPartStart);

                    %for the last time, get switching points and class durations
                    switch123 = getSwitchPointClassDuration(pointClassMSS);
                    switch123 = switch123(:,1);

                    %assign subpart classification
                    partClassMSS = [switch123(1:end-1) switch123(2:end)-1 pointClassMSS(switch123(1:end-1))];

                    %shift time to start part classification at part start time
                    partClassMSS(:,1:2) = partClassMSS(:,1:2) + trackPartStart - 1;

                    %get number of subparts
                    numSubparts = size(partClassMSS,1);

                    %go over subparts and get their classification
                    classMSSTmp = NaN(numSubparts,1);
                    for iSubpart = 1 : numSubparts

                        %get subpart's start and end points
                        startPoint = 8*(partClassMSS(iSubpart,1)-1)+1;
                        endPoint = 8*partClassMSS(iSubpart,2);

                        %call MSS analysis code
                        classMSSTmp(iSubpart) = trackMSSAnalysis(tracks(iTrack,startPoint:endPoint),...
                            probDim,momentOrders,alphaValues);

                    end

                    %update subpart classification
                    classMSSTmp(isnan(classMSSTmp)) = 0;
                    partClassMSS(:,3) = classMSSTmp;

                    %merge subparts that now have the same classification
                    for iSubpart = 1 : numSubparts-1
                        iSubpartPlus1 = iSubpart + 1;
                        while ( (iSubpartPlus1 <= numSubparts) && ...
                                ( (classMSSTmp(iSubpart) == classMSSTmp(iSubpartPlus1)) || ...
                                (isnan(classMSSTmp(iSubpart)) && isnan(classMSSTmp(iSubpartPlus1))) ) )
                            partClassMSS(iSubpart,2) = partClassMSS(iSubpartPlus1,2);
                            partClassMSS(iSubpartPlus1,1) = partClassMSS(iSubpart,1);
                            iSubpartPlus1 = iSubpartPlus1 + 1;
                        end
                    end
                    [dummy,uniqueParts] = unique(partClassMSS(:,1));
                    partClassMSS = partClassMSS(uniqueParts,:);
                    numSubparts = size(partClassMSS,1);

                    %now go over subparts and get final diffusion characteristics
                    for iSubpart = 1 : numSubparts

                        %get subpart's start and end points
                        startPoint = 8*(partClassMSS(iSubpart,1)-1)+1;
                        endPoint = 8*partClassMSS(iSubpart,2);

                        %call MSS analysis code
                        [classMSSTmp,mssSlopeTmp,genDiffCoefTmp,scalingPowerTmp,normDiffCoefTmp] = ...
                            trackMSSAnalysis(tracks(iTrack,startPoint:endPoint),...
                            probDim,momentOrders,alphaValues);

                        %estimate confinement radius (and center for plotting)
                        if classMSSTmp == 1
                            [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,startPoint:endPoint),probDim,confRadMin);
                        else
                            confRadTmp = NaN;
                            centerTmp = NaN(1,probDim);
                        end

                        %store this subpart's information
                        partClassMSS(iSubpart,3:20+probDim) = [classMSSTmp ...
                            mssSlopeTmp genDiffCoefTmp scalingPowerTmp ...
                            normDiffCoefTmp confRadTmp centerTmp];

                    end

                else

                    %if there are no classified points, then keep whole part as
                    %unclassified
                    partClassMSS = [trackPartStart trackPartStart+trackPartLength(iPart)-1 NaN(1,18+probDim)];

                end %(if ~isempty(pointClassMSS~=0) ... else ...)
                
            else %if number of rolling windows is not larger than minimum duration, classify track part as a whole

                %get part's start and end points
                partClassMSS = [];
                partClassMSS(1,1:2) = [trackPartStart trackPartStart+trackPartLength(iPart)-1];
                startPoint = 8*(partClassMSS(1,1)-1)+1;
                endPoint = 8*partClassMSS(1,2);
                
                %call MSS analysis code
                [classMSSTmp,mssSlopeTmp,genDiffCoefTmp,scalingPowerTmp,normDiffCoefTmp] = ...
                    trackMSSAnalysis(tracks(iTrack,startPoint:endPoint),...
                    probDim,momentOrders,alphaValues);

                %estimate confinement radius (and center for plotting)
                %estimate confinement radius (and center for plotting)
                if classMSSTmp == 1
                    [confRadTmp,centerTmp] = estimConfRad(tracks(iTrack,startPoint:endPoint),probDim,confRadMin);
                else
                    confRadTmp = NaN;
                    centerTmp = NaN(1,probDim);
                end

                %store this subpart's information
                partClassMSS(1,3:20+probDim) = [classMSSTmp mssSlopeTmp ...
                    genDiffCoefTmp scalingPowerTmp normDiffCoefTmp ...
                    confRadTmp centerTmp];

            end %(if numRollWindows > minDuration(1))

            %get number of subparts for this track segment part
            numSubparts = size(partClassMSS,1);

            %insert the new subparts in their proper location in the track
            %segment classification
            trackClassMSS = [trackClassMSS(1:oldPart2newPartMap(iPart)-1,:); ...
                partClassMSS; trackClassMSS(oldPart2newPartMap(iPart)+1:end,:)];

            %also divide the expanded asymmetry classification into the
            %same subparts
            asymClass = trackClassAsymExp(oldPart2newPartMap(iPart),3);
            asymInsertParts = [partClassMSS(:,1:2) repmat(asymClass,numSubparts,1)];
            trackClassAsymExp = [trackClassAsymExp(1:oldPart2newPartMap(iPart)-1,:); ...
                asymInsertParts; trackClassAsymExp(oldPart2newPartMap(iPart)+1:end,:)];

            %update the mapping from old parts to new parts
            oldPart2newPartMap(iPart+1:end) = oldPart2newPartMap(iPart+1:end) + numSubparts - 1;

        end %(for iPart = 1 : numParts2analyze)

    end %(if ~isempty(tracksParts2analyze))

    %% Store track segment information
    trackSegmentClassRes(iTrack).asymmetry = trackClassAsym;
    trackSegmentClassRes(iTrack).momentScalingSpectrum = trackClassMSS;
    trackSegmentClassRes(iTrack).asymmetryAfterMSS = trackClassAsymExp;

end %(for iTrack = indx4analysis')

%% Store trivial nonclassification information for tracks that are not classifiable
for iTrack = indxNot4analysis'
    trackSELCurrent = trackSEL(iTrack,:);
    trackSegmentClassRes(iTrack).asymmetry(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).momentScalingSpectrum(1:2) = trackSELCurrent(1:2);
    trackSegmentClassRes(iTrack).asymmetryAfterMSS(1:2) = trackSELCurrent(1:2);
end

%% save results in output structure

%reserve memory
segmentClass = struct('asymmetry',[],'momentScalingSpectrum',[],'asymmetryAfterMSS',[]);
transDiffAnalysisRes = repmat(struct('segmentClass',segmentClass),numInputTracks,1);

%go over all input tracks
for iTrack = 1 : numInputTracks

    %go over the segments of each track
    for iSegment = 1 : numSegments(iTrack)

        %store the segment's classification results
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetry = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetry;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).momentScalingSpectrum = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).momentScalingSpectrum;
        transDiffAnalysisRes(iTrack).segmentClass(iSegment,1).asymmetryAfterMSS = ...
            trackSegmentClassRes(compTrackStartRow(iTrack)+iSegment-1).asymmetryAfterMSS;

    end %(for iSegment = 1 : numSegments(iTrack))

end %(for iTrack = 1 : numInputTracks)

%% plotting

%plot results if requested
if plotRes
    plotTracksTransDiffAnalysis2D(tracksInput,transDiffAnalysisRes,[],1);
end


%% Subfunctions

function pointClassMSS = reclassUnclassPoints(pointClassMSS,mssSlope,...
    tracks,halfWindowMSS,probDim,momentOrders,alphaMSS,expandWindow,...
    trackPartStart)

%get the switching points and class durations
[switch123,classDuration] = getSwitchPointClassDuration(pointClassMSS);

%shift switching times by track start time
switch123Mod(:,1) = switch123(:,1) + trackPartStart + halfWindowMSS*expandWindow - 1;

%find all subparts that are unclassified
badClass = find(classDuration(:,2) == 0);

if ~isempty(badClass) && size(classDuration,1) > 1

    %if the first subpart is at the track start, lump it with the subpart
    %after and reclassify the merged subparts
    %if the merged subparts get the same classification as the second
    %subpart alone, keep them merged
    %otherwise, keep unclassified subpart as unclassified
    if badClass(1) == 1
        lumpedPartsStart = switch123Mod(1,1) - halfWindowMSS*expandWindow;
        lumpedPartsEnd = switch123Mod(3,1) - 1 + halfWindowMSS*expandWindow;
        lumpedPartsClass = trackMSSAnalysis(tracks(...
            8*(lumpedPartsStart-1)+1:8*lumpedPartsEnd),...
            probDim,momentOrders,alphaMSS);
        if lumpedPartsClass == classDuration(2,2)
            pointClassMSS(switch123(1,1):switch123(2,1)-1) = ...
                classDuration(2,2);
        end
        badClass = badClass(2:end);
    end
    
    if ~isempty(badClass)
        
        %if the last subpart is at the track end, lump it with the subpart
        %before and reclassify the merged subparts
        %if the merged subparts get the same classification as the
        %second-to-last subpart alone, keep them merged
        %otherwise, keep unclassified subpart as unclassified
        if badClass(end) == size(classDuration,1)
            lumpedPartsStart = switch123Mod(end-2,1) - halfWindowMSS*expandWindow;
            lumpedPartsEnd = switch123Mod(end,1) - 1 + halfWindowMSS*expandWindow;
            lumpedPartsClass = trackMSSAnalysis(tracks(...
                8*(lumpedPartsStart-1)+1:8*lumpedPartsEnd),...
                probDim,momentOrders,alphaMSS);
            if lumpedPartsClass == classDuration(end-1,2)
                pointClassMSS(switch123(end-1,1):switch123(end,1)-1) = ...
                    classDuration(end-1,2);
            end
            badClass = badClass(1:end-1);
        end
        
        %if there are subparts in the middle ...
        if ~isempty(badClass)
            
            %go over each subpart
            for iClass = badClass'

                %if the subparts before and after have the same
                %classification, lump all three subparts together and
                %reclassify the merged subparts
                %if the merged subparts get the same classification as the
                %separate subparts, then keep them merged
                %otherwise, keep unclassified subpart as unclassified
                if classDuration(iClass-1,2) == classDuration(iClass+1,2)
                    lumpedPartsStart = switch123Mod(iClass-1,1) - ...
                        halfWindowMSS*expandWindow;
                    lumpedPartsEnd = switch123Mod(iClass+2,1) - 1 + ...
                        halfWindowMSS*expandWindow;
                    lumpedPartsClass = trackMSSAnalysis(tracks(...
                        8*(lumpedPartsStart-1)+1:8*lumpedPartsEnd),...
                        probDim,momentOrders,alphaMSS);
                    if lumpedPartsClass == classDuration(iClass-1,2)
                        pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = ...
                            classDuration(iClass-1,2);
                    end
                    
                else %if the subparts before and after have different classifications,
                    %lump the unclassified subpart first with the one
                    %before and then with the one after it, reclassify each
                    %combination, and then reclassify the unclassified
                    %subpart according to which reclassification is more
                    %consistent with previous classification

                    %with subpart before
                    lumpedPartsStart = switch123Mod(iClass-1,1) - ...
                        halfWindowMSS*expandWindow;
                    lumpedPartsEnd = switch123Mod(iClass+1,1) - 1 + ...
                        halfWindowMSS*expandWindow;
                    [lumpedPartsClass1,mssSlopeLumped1] = trackMSSAnalysis(...
                        tracks(8*(lumpedPartsStart-1)+1:8*lumpedPartsEnd),...
                        probDim,momentOrders,alphaMSS);
                    
                    %with subpart after
                    lumpedPartsStart = switch123Mod(iClass,1) - ...
                        halfWindowMSS*expandWindow;
                    lumpedPartsEnd = switch123Mod(iClass+2,1) - 1 + ...
                        halfWindowMSS*expandWindow;
                    [lumpedPartsClass2,mssSlopeLumped2] = trackMSSAnalysis(...
                        tracks(8*(lumpedPartsStart-1)+1:8*lumpedPartsEnd),...
                        probDim,momentOrders,alphaMSS);
                    
                    %compare the two classifications with previous
                    %classifications
                    if lumpedPartsClass1 == classDuration(iClass-1,2)
                        if lumpedPartsClass2 == classDuration(iClass+1,2)
                            %agree with both - most complicated case
                            mssSlopeClassMinus1 = nanmean(...
                                mssSlope(switch123(iClass-1,1):...
                                switch123(iClass,1)-1));
                            mssSlopeClassPlus1 = nanmean(...
                                mssSlope(switch123(iClass+1,1):...
                                switch123(iClass+2,1)-1));
                            slopeDiff1 = abs(mssSlopeLumped1-...
                                mssSlopeClassMinus1);
                            slopeDiff2 = abs(mssSlopeLumped2-...
                                mssSlopeClassPlus1);
                            if slopeDiff1 <= slopeDiff2
                                pointClassMSS(switch123(iClass,...
                                    1):switch123(iClass+1,1)-1)...
                                    = classDuration(iClass-1,2);
                            else
                                pointClassMSS(switch123(iClass,...
                                    1):switch123(iClass+1,1)-1)...
                                    = classDuration(iClass+1,2);
                            end
                        else
                            %agree only with subpart before
                            pointClassMSS(switch123(iClass,...
                                1):switch123(iClass+1,1)-1)...
                                = classDuration(iClass-1,2);
                        end
                    else
                        if lumpedPartsClass2 == classDuration(iClass+1,2)
                            %agree only with subpart after
                            pointClassMSS(switch123(iClass,...
                                1):switch123(iClass+1,1)-1)...
                                = classDuration(iClass+1,2);
                        end
                    end %(if lumpedPartsClass1 == classDuration(iClass-1,2))
                    
                end %(if classDuration(iClass-1,2) == classDuration(iClass+1,2))
                
            end %(for iClass = badClass')
            
        end %(if ~isempty(badClass))
        
    end %(if ~isempty(badClass))
    
end %(if ~isempty(badClass))



function [pointClassMSS,mssSlope,halfWindowStartRemoved,halfWindowEndRemoved] ...
    = reclassConfinedDirectedPoints(pointClassMSS,mssSlope,minDuration,halfWindowMSS)

%get the switching points and class durations
[switch123,classDuration] = getSwitchPointClassDuration(pointClassMSS);

%find all subparts originally classified as confined or directed but for 
%less than their respective minDuration frames, and convert them into pure Brownian
badClass = find( (classDuration(:,1)<minDuration(1) & classDuration(:,2)==1) ...
    | (classDuration(:,1)<minDuration(2) & classDuration(:,2)==3) );
for iClass = badClass'
    pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = 2;
end

%add to pointClassMSS and mssSlope the first and last halfWindowMSS points
pointClassMSS = [zeros(halfWindowMSS,1); pointClassMSS; zeros(halfWindowMSS,1)];
mssSlope = [NaN(halfWindowMSS,1); mssSlope; NaN(halfWindowMSS,1)];

%get points still classified as confined or directed
pointClassMSSTmp = pointClassMSS;
indxGood = find(pointClassMSS==1 | pointClassMSS==3);

%go over these points and expand their classification to the points around 
%them and constituting their window
for iPoint = indxGood'

    %get point classification, thus determine the opposite classification
    currentClass = pointClassMSSTmp(iPoint);
    oppositeClass = setdiff([1 3],currentClass);

    %find the start and end points of this point's window
    pointMin = iPoint - halfWindowMSS;
    pointMax = iPoint + halfWindowMSS;

    %find the closest points of opposite classification on either side
    oppClassMin = find(pointClassMSSTmp(1:iPoint-1)==oppositeClass,1,'last');
    oppClassMax = find(pointClassMSSTmp(iPoint+1:end)==oppositeClass,1,'first') ...
        + iPoint;

    %get the windows of those points
    pointMinMax = oppClassMin + halfWindowMSS;
    pointMaxMin = oppClassMax - halfWindowMSS;

    %if there is any overlap between these windows and the current point's
    %window, modify the current point's window
    if pointMinMax > pointMin
        pointMin = ceil(mean([pointMinMax pointMin]));
    end
    if pointMaxMin < pointMax
        pointMax = floor(mean([pointMaxMin pointMax]));
    end

    %now expand the classification
    pointClassMSS(pointMin:pointMax) = pointClassMSSTmp(iPoint);

end

%remove the first and last halfWindowMSS points if they remained unclassified
halfWindowStartRemoved = 0;
if pointClassMSS(1) == 0
    pointClassMSS = pointClassMSS(halfWindowMSS+1:end);
    mssSlope = mssSlope(halfWindowMSS+1:end);
    halfWindowStartRemoved = 1;
end
halfWindowEndRemoved = 0;
if pointClassMSS(end) == 0
    pointClassMSS = pointClassMSS(1:end-halfWindowMSS);
    mssSlope = mssSlope(1:end-halfWindowMSS);
    halfWindowEndRemoved = 1;
end



function pointClassMSS = reclassBrownianPoints(pointClassMSS,mssSlope,windowBrown)

%get the switching points and class durations
[switch123,classDuration] = getSwitchPointClassDuration(pointClassMSS);

%find all subparts originally classified as Brownian, but for
%less than windowBrown frames
badClass = find(classDuration(:,1) < windowBrown & classDuration(:,2) == 2);

if ~isempty(badClass) && size(classDuration,1) > 1
    
    %if the first subpart is at the track start, give it the classification
    %of the subpart just after it
    if badClass(1) == 1
        pointClassMSS(switch123(1,1):switch123(2,1)-1) = classDuration(2,2);
        badClass = badClass(2:end);
    end
    
    if ~isempty(badClass)
        
        %if the last subpart is at the track end, give it the classification 
        %of the subpart just before it
        if badClass(end) == size(classDuration,1)
            pointClassMSS(switch123(end-1,1):switch123(end,1)-1) = classDuration(end-1,2);
            badClass = badClass(1:end-1);
        end
        
        %if there are subparts in the middle ...
        if ~isempty(badClass)
            
            %go over each subpart
            for iClass = badClass'
                
                %if the subparts before and after have the same
                %classification, give the previously-Brownian subpart that
                %classification
                if classDuration(iClass-1,2) == classDuration(iClass+1,2)
                    pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = classDuration(iClass-1,2);

                else %if the subparts before and after have different classifications,
                    %calculate the average slope of the previouly-Brownian
                    %subpart, and assign it a temporary classification
                    averageSlope = nanmean(mssSlope(switch123(iClass,1):switch123(iClass+1,1)-1));
                    pointClassMSSBB = 2*(averageSlope > 0.5) + 1;
                    
                    %if the temporary classification is the same as either
                    %adjacent subpart classification, retain the classification
                    if any(classDuration([iClass-1;iClass+1],2) == pointClassMSSBB)
                        pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = pointClassMSSBB;
                    else
                        pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = 0;
                    end
                    
                end %(if classDuration(iClass-1,2) == classDuration(iClass+1,2))
                
            end %(for iClass = badClass')
            
        end %(if ~isempty(badClass))
        
    end %(if ~isempty(badClass))
    
end %(if ~isempty(badClass))



function [switch123,classDuration] = getSwitchPointClassDuration(pointClassMSS)

%take the pointClassMSS difference: 0 means staying in the same state,
%while nonzero values indicate motion type switches
%add one point to the beginning and the end, to get end points
pointClassMSSFirst = setdiff([1 2 3],pointClassMSS(1));
pointClassMSSFirst = pointClassMSSFirst(1);
pointClassMSSLast = setdiff([1 2 3],pointClassMSS(end));
pointClassMSSLast = pointClassMSSLast(1);
pointClassMSSDiff = diff([pointClassMSSFirst; pointClassMSS; pointClassMSSLast]);
switch123 = find(pointClassMSSDiff ~= 0);

%make an array of classes and their starting points note that the very last
%class is not real
switch123 = [switch123 [pointClassMSS(switch123(1:end-1)); pointClassMSSLast]];

%get the duration of each class (and with this ignore the imaginary last
%class)
classDuration = [diff(switch123(:,1)) switch123(1:end-1,2)];


function [confRadTmp,centerTmp] = estimConfRad(tracks,probDim,confRadMin)

%get subpart's coordinates
xCoord = tracks(1:8:end)';
yCoord = tracks(2:8:end)';
zCoord = tracks(3:8:end)';
xyzCoord = [xCoord yCoord zCoord];

%find the eigenvalues and eigenvectors of the variance-covariance
%matrix of this track's positions
eigenVal = eig(nancov(xyzCoord(:,1:probDim)));

%calculate the track's confinement radius
if confRadMin
    confRadTmp = sqrt( min(eigenVal) * (probDim + 2) );
else
    confRadTmp = sqrt( mean(eigenVal) * (probDim + 2) );
end

%calculate the track's center
centerTmp = nanmean(xyzCoord(:,1:probDim));
