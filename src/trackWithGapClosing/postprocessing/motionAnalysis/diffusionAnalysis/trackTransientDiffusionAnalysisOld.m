function [transDiffAnalysisRes,errFlag] = trackTransientDiffusionAnalysisOld(tracks,...
    extractType,probDim,checkAsym,alphaValues,minDuration,plotRes,confRadMin)
%TRACKTRANSIENTDIFFUSIONANALYSIS performs a rolling window diffusion analysis, checking first for asymmetry
%
%SYNOPSIS [transDiffAnalysisRes,errFlag] = trackTransientDiffusionAnalysis1(tracks,...
%    extractType,probDim,checkAsym,alphaValues,minDuration,plotRes,confRadMin)
%
%INPUT  tracks      : -- EITHER --
%                     Output of trackWithGapClosing (matrix),
%                     -- OR --
%                     Output of trackCloseGapsKalman (structure, possibly
%                     with merges/splits.
%       extractType : 1 - Analyze every track segment separately.
%                     2 - Extract from each compound track the longest
%                         trajectory to use in analysis - NOT IMPLEMENTED
%                         YET.
%                     Variable irrelevant if tracks are input as a matrix.
%                     Optional. Default: 1.
%       probDim     : Problem dimensionality.
%                     Optional. Default: 2.
%       checkAsym   : 1 to check for asymmetric tracks and to analyze their
%                     diffusion after dimensionality reduction, 0
%                     otherwise.
%                     Optional. Default: 0.
%       alphaValues : Row vector with 2 entries. First entry is the
%                     alpha-value for MSS analysis (can take the values
%                     0.2, 0.1, 0.05 and 0.01). Second entry is the
%                     alpha-value for asymmetry determination (can take the
%                     values 0.2, 0.1 and 0.05).
%                     Optional. Default: [0.1 0.1]. If only one value is
%                     entered, it is taken as the alpha-value for MSS
%                     analysis.
%       minDuration : Row vector with 3 entries, with minimum durations for
%                     confined, directed, and asymmetric trajectory
%                     segments. Last entry not needed if checkAsym = 0.
%                     Must be greater than or equal to [20 20 5].
%                     Optional. Default: [20 20 5].
%       plotRes     : 1 to plot results, 0 otherwise.
%                     Optional. Default: 0.
%                     Results can be plotted only if problem is 2D.
%                     color-coding:
%                     *red: linear via asymmetry.
%                     *blue: random/unclassified & confined diffusion via MSS.
%                     *cyan: random/unclassified & normal diffusion via MSS.
%                     *magenta: random/unclassified & super diffusion via MSS.
%                     *purple: random via asymmetry & unclassifiable via MSS.
%                     *black: unclassified at all.
%       confRadMin  : 1 to calculate the confinement radius of confined
%                     particles using the minimum positional standard
%                     deviation, 0 to calculate it using the mean
%                     positional standard deviation.
%                     Optional. Default: 0.
%
%OUTPUT transDiffAnalysisRes : ...
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
    disp('--trackDiffusionAnalysis1: Please input at least the tracks to be analyzed!');
    errFlag = 1;
    return
end

if nargin < 2 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--trackDiffusionAnalysis1: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(checkAsym)
    checkAsym = 0;
end

if nargin < 5 || isempty(alphaValues)
    alphaValues = [0.1 0.1];
elseif length(alphaValues) == 1
    alphaValues = [alphaValues 0.1];
end
alphaMSS = alphaValues(1);
alphaAsym = alphaValues(2);

if nargin < 6 || isempty(minDuration)
    minDuration = [20 20 5];
else
    if length(minDuration) == 2
        minDuration = [minDuration 5];
    end
    if minDuration(1) < 20
        minDuration(1) = 20;
    end
    if minDuration(2) < 20
        minDuration(2) = 20;
    end
    if minDuration(3) < 5
        minDuration(3) = 5;
    end
end

if nargin < 7 || isempty(plotRes)
    plotRes = 0;
elseif plotRes == 1 && probDim ~= 2
    disp('--trackDiffusionAnalysis1: Cannot plot tracks if problem is not 2D!');
    plotRes = 0;
end

if nargin < 8 || isempty(confRadMin)
    confRadMin = 0;
end

if errFlag
    disp('--trackDiffusionAnalysis1: Please fix input variables');
    return
end

%define window sizes
windowAsym = 5;
windowMSS = 21;

%define duration of a Brownian segment the leads to lumping it with the
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
    indx4analysis = find(trackSEL(:,3) >= windowMSS);
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

        %get number of asymmetry analysis rolling windows for this track
        %segment
        numRollWindows = trackSELCurrent(3) - windowAsym + 1;

        %initialize point classification vector
        pointClassAsym = NaN(trackSELCurrent(3),1);

        %get the particle positions along the track segment
        coordX = tracks(iTrack,1:8:end)';
        coordY = tracks(iTrack,2:8:end)';
        coordZ = tracks(iTrack,3:8:end)';
        coordXYZ = [coordX coordY coordZ];
        coordXYZ = coordXYZ(trackSELCurrent(1):trackSELCurrent(2),:);

        %go over all windows
        for iWindow = 1 : numRollWindows

            %determine whether the track segment is sufficiently asymmetric
            %in this window
            coord2classify = coordXYZ(iWindow:iWindow+windowAsym-1,1:probDim);
            if any(~isnan(coord2classify))
                [asymParamT,asymFlag] = asymDeterm2D3D(coord2classify,alphaAsym);
            else
                asymFlag = NaN;
            end

            %classify window middle point as ...
            %1 = linear, if the asymmetry parameter is larger than the threshold
            %0 = not linear, if the asymmetry parameter is smaller than the
            %threshold
            %otherwise, keep track segment classification as undetermined
            pointClassAsym(iWindow+(windowAsym-1)/2) = asymFlag;

        end

        %for any point classified as asymmetric, expand its classification
        %to the points surrounding it and corresponding to its window
        indxGood = find(pointClassAsym==1);
        for iPoint = indxGood'
            pointClassAsym(iPoint-(windowAsym-1)/2:iPoint+(windowAsym-1)/2) = 1;
        end

        %expand the classification to the first and last few points of a
        %track segment
        pointClassAsym(1:(windowAsym-1)/2) = pointClassAsym((windowAsym+1)/2);
        pointClassAsym(end-(windowAsym-3)/2:end) = pointClassAsym(end-(windowAsym-1)/2);

        %give unclassified points the classification of the point before them
        unclassPoints = find(isnan(pointClassAsym));
        if ~isempty(unclassPoints)
            if unclassPoints(1) == 1
                unclassPoints = unclassPoints(2:end);
            end
            for iPoint = unclassPoints'
                pointClassAsym(iPoint) = pointClassAsym(iPoint-1);
            end
        end

        %if this track segment got some classification ...
        indxGood = find(~isnan(pointClassAsym));
        if ~isempty(indxGood)

            %if there are unclassified points in the beginning, give them the
            %classification of the first classified point
            if isnan(pointClassAsym(1))
                pointClassAsym(1:indxGood(1)-1) = pointClassAsym(indxGood(1));
            end

            %take the pointClassAsym difference: 0 means staying in the same
            %state, 1 means switching from symmetric to asymmetric, and -1
            %means switching from asymmetric to symmetric
            %add one point to the beginning and the end, to get end points
            pointClassAsymDiff = diff([1-pointClassAsym(1); pointClassAsym; 1-pointClassAsym(end)]);
            switch0110 = find(pointClassAsymDiff ~= 0);

            %make an array of classes and their starting points
            %note that the very last class is not real
            switch0110 = [switch0110 [pointClassAsym(switch0110(1:end-1)); 1-pointClassAsym(end)]]; %#ok

            %get the duration of each class
            classDuration = [diff(switch0110(:,1)) switch0110(1:end-1,2)];

            %find all track segment parts originally classified as asymmetric but
            %for less than minDuration(3) frames, and convert them into symmetric
            badClass = find(classDuration(:,1) < minDuration(3) & classDuration(:,2)==1);
            for iClass = badClass'
                pointClassAsym(switch0110(iClass,1):switch0110(iClass+1,1)-1) = 0;
            end

            %take the pointClassAsym difference again and get class durations
            %again ...
            pointClassAsymDiff = diff([1-pointClassAsym(1); pointClassAsym; 1-pointClassAsym(end)]);
            switch0110 = find(pointClassAsymDiff ~= 0);
            switch0110 = [switch0110 [pointClassAsym(switch0110(1:end-1)); 1-pointClassAsym(end)]]; %#ok
            classDuration = [diff(switch0110(:,1)) switch0110(1:end-1,2)];

            %find all track segment parts originally classified as
            %symmetric, but for less than windowAsym frames, and convert
            %them into asymmetric
            badClass = find(classDuration(:,1) < windowAsym & classDuration(:,2)==0);
            for iClass = badClass'
                pointClassAsym(switch0110(iClass,1):switch0110(iClass+1)-1,1) = 1;
            end

            %for the last time, find switching points
            %this is the final asymmetry classification of this track
            %segment
            pointClassAsymDiff = diff([1-pointClassAsym(1); pointClassAsym; 1-pointClassAsym(end)]);
            switch0110 = find(pointClassAsymDiff ~= 0);
            trackClassAsym = [switch0110(1:end-1) switch0110(2:end)-1 pointClassAsym(switch0110(1:end-1))];

            %shift time to start track segment classification at track
            %segment start time
            trackClassAsym(:,1:2) = trackClassAsym(:,1:2) + trackSELCurrent(1) - 1;

        else

            %if there are no classified points, then keep whole track
            %segment as unclassified
            trackClassAsym = [trackSELCurrent(1:2) NaN];

        end %(if ~isempty(indxGood) ... else ...)

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
    %that are longer than windowMSS
    trackParts2analyze = find(trackClassMSS(:,3) ~= -1 & trackPartLength >= windowMSS);

    %if there are parts to analyze ...
    if ~isempty(trackParts2analyze)

        %go over these parts ...
        for iPart = trackParts2analyze'

            %get starting point of this part
            trackPartStart = trackClassMSS(oldPart2newPartMap(iPart),1);

            %get number of MSS analysis rolling windows in this part
            numRollWindows = trackPartLength(iPart) - windowMSS + 1;

            %initialize point classification vector
            pointClassMSS = NaN(trackPartLength(iPart),1);
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
                [pointClassMSS(iWindow+(windowMSS-1)/2),mssSlope(iWindow+(windowMSS-1)/2)] = ...
                    trackMSSAnalysis(tracks(iTrack,8*(startPoint-1)+1:8*endPoint),probDim,momentOrders,alphaMSS);

            end %(for iWindow = 1 : numRollWindows)

            %for any point classified as confined or directed, expand its
            %classification to the points surrounding it and corresponding
            %to its window ...

            %             %first get points classified as confined or directed
            %             pointClassMSSTmp = pointClassMSS;
            %             mssSlopeTmp = mssSlope;
            %             indxGood = find(pointClassMSS==1 | pointClassMSS==3);
            %
            %             %go over these points
            %             for iPoint = indxGood'
            %
            %                 %get point classification, thus determine the opposite
            %                 %classification
            %                 currentClass = pointClassMSSTmp(iPoint);
            %                 oppositeClass = setdiff([1 3],currentClass);
            %
            %                 %find the start and end points of this point's window
            %                 pointMin = iPoint - (windowMSS-1)/2;
            %                 pointMax = iPoint + (windowMSS-1)/2;
            %
            %                 %find the closest points of opposite classification on
            %                 %either side
            %                 oppClassMin = find(pointClassMSSTmp(1:iPoint-1)==oppositeClass,1,'last');
            %                 oppClassMax = find(pointClassMSSTmp(iPoint+1:end)==oppositeClass,1,'first') ...
            %                     + iPoint;
            %
            %                 %get the windows of those points
            %                 pointMinMax = oppClassMin + (windowMSS-1)/2;
            %                 pointMaxMin = oppClassMax - (windowMSS-1)/2;
            %
            %                 %if there is any overlap between these windows and the
            %                 %current point's window, modify the current point's window
            %                 if pointMinMax > pointMin
            %                     pointMin = ceil(mean([pointMinMax pointMin]));
            %                 end
            %                 if pointMaxMin < pointMax
            %                     pointMax = floor(mean([pointMaxMin pointMax]));
            %                 end
            %
            %                 %now expand classification
            %                 pointClassMSS(pointMin:pointMax) = pointClassMSSTmp(iPoint);
            %                 mssSlope(pointMin:pointMax) = mssSlopeTmp(iPoint);
            %             end

            %             %expand the classification to the first and last few points of a
            %             %track segment
            %             pointClassMSS(1:(windowMSS-1)/2) = pointClassMSS((windowMSS+1)/2);
            %             pointClassMSS(end-(windowMSS-3)/2:end) = pointClassMSS(end-(windowMSS-1)/2);
            %             mssSlope(1:(windowMSS-1)/2) = mssSlope((windowMSS+1)/2);
            %             mssSlope(end-(windowMSS-3)/2:end) = mssSlope(end-(windowMSS-1)/2);

            %give unclassified points the classification of the point before them
            unclassPoints = find(isnan(pointClassMSS));
            if ~isempty(unclassPoints)
                if unclassPoints(1) == 1
                    unclassPoints = unclassPoints(2:end);
                end
                for iPoint = unclassPoints'
                    pointClassMSS(iPoint) = pointClassMSS(iPoint-1);
                    mssSlope(iPoint) = mssSlope(iPoint-1);
                end
            end

            %if this part got some classification ...
            indxGood = find(~isnan(pointClassMSS));
            if ~isempty(indxGood)

                %if there are unclassified points in the beginning, give them the
                %classification of the first classified point
                if isnan(pointClassMSS(1))
                    pointClassMSS(1:indxGood(1)-1) = pointClassMSS(indxGood(1));
                    mssSlope(1:indxGood(1)-1) = mssSlope(indxGood(1));
                end

                %take the pointClassMSS difference: 0 means staying in the same
                %state, 1 & 2 mean switching to a state of higher mobility,
                %-1 & -2 mean switching to a state of lower mobility
                %add one point to the beginning and the end, to get end points
                pointClassMSSFirst = setdiff([1 2 3],pointClassMSS(1));
                pointClassMSSFirst = pointClassMSSFirst(1);
                pointClassMSSLast = setdiff([1 2 3],pointClassMSS(end));
                pointClassMSSLast = pointClassMSSLast(1);
                pointClassMSSDiff = diff([pointClassMSSFirst; pointClassMSS; pointClassMSSLast]);
                switch123 = find(pointClassMSSDiff ~= 0);

                %make an array of classes and their starting points
                %note that the very last class is not real
                switch123 = [switch123 [pointClassMSS(switch123(1:end-1)); pointClassMSSLast]]; %#ok

                %get the duration of each class
                classDuration = [diff(switch123(:,1)) switch123(1:end-1,2)];

                %find all subparts originally classified as confined or
                %directed but for less than their corresponding minDuration
                %frames, and convert them into pure Brownian
                badClass = find( (classDuration(:,1)<minDuration(1) & classDuration(:,2)==1) ...
                    | (classDuration(:,1)<minDuration(2) & classDuration(:,2)==3) );
                for iClass = badClass'
                    pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = 2;
                end

                %take the pointClassMSS difference again and get class durations
                %again ...
                pointClassMSSFirst = setdiff([1 2 3],pointClassMSS(1));
                pointClassMSSFirst = pointClassMSSFirst(1);
                pointClassMSSLast = setdiff([1 2 3],pointClassMSS(end));
                pointClassMSSLast = pointClassMSSLast(1);
                pointClassMSSDiff = diff([pointClassMSSFirst; pointClassMSS; pointClassMSSLast]);
                switch123 = find(pointClassMSSDiff ~= 0);
                switch123 = [switch123 [pointClassMSS(switch123(1:end-1)); pointClassMSSLast]]; %#ok
                classDuration = [diff(switch123(:,1)) switch123(1:end-1,2)];

                %find all subparts originally classified as Brownian, but for
                %less than windowBrown frames
                %if the subparts before and after such a subpart have the
                %same classification, give that classification to the
                %previously-Brownian subpart
                %if the subparts before and after have different
                %classifications, then keep segment as separate
                badClass = find(classDuration(:,1) < windowBrown & classDuration(:,2) == 2);
                if ~isempty(badClass)
                    if badClass(1) == 1
                        pointClassMSS(switch123(1,1):switch123(2,1)-1) = classDuration(2,2);
                        badClass = badClass(2:end);
                    end
                    if ~isempty(badClass)
                        if badClass(end) == size(classDuration,1);
                            pointClassMSS(switch123(end-1,1):switch123(end,1)-1) = classDuration(end-1,2);
                            badClass = badClass(1:end-1);
                        end
                        if ~isempty(badClass)
                            for iClass = badClass'
                                if classDuration(iClass-1,2) == classDuration(iClass+1,2)
                                    pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = classDuration(iClass-1,2);
                                else
                                    averageSlope = mean(mssSlope(switch123(iClass,1):switch123(iClass+1,1)-1));
                                    pointClassMSS(switch123(iClass,1):switch123(iClass+1,1)-1) = 2*(averageSlope > 0.5) + 1;
                                end
                            end
                        end
                    end
                end

                %for the last time, find switching points
                %this is the final MSS classification of this track segment
                %part
                pointClassMSSFirst = setdiff([1 2 3],pointClassMSS(1));
                pointClassMSSFirst = pointClassMSSFirst(1);
                pointClassMSSLast = setdiff([1 2 3],pointClassMSS(end));
                pointClassMSSLast = pointClassMSSLast(1);
                pointClassMSSDiff = diff([pointClassMSSFirst; pointClassMSS; pointClassMSSLast]);
                switch123 = find(pointClassMSSDiff ~= 0);
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
                        probDim,momentOrders,alphaMSS);

                end

                %update subpart classification
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
                        probDim,momentOrders,alphaMSS);

                    %estimate confinement radius (and center for plotting)
                    if classMSSTmp == 1

                        %get subpart's coordinates
                        xCoord = (tracks(iTrack,startPoint:8:endPoint))';
                        yCoord = (tracks(iTrack,startPoint+1:8:endPoint))';
                        zCoord = (tracks(iTrack,startPoint+2:8:endPoint))';
                        xyzCoord = [xCoord yCoord zCoord];

                        %find the eignevalues and eigenvectors of the variance-covariance
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

                    else
                        confRadTmp = NaN;
                        centerTmp = NaN(1,probDim);
                    end

                    %store this subpart's information
                    partClassMSS(iSubpart,3:(20+probDim)) = [classMSSTmp ...
                        mssSlopeTmp genDiffCoefTmp scalingPowerTmp ...
                        normDiffCoefTmp confRadTmp centerTmp];

                end

            else

                %if there are no classified points, then keep whole part as
                %unclassified
                partClassMSS = [trackPartStart trackPartStart+trackPartLength(iPart)-1 NaN(1,18+probDim)];

            end %(if ~isempty(indxGood) ... else ...)

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

%% ~~~ the end ~~~
