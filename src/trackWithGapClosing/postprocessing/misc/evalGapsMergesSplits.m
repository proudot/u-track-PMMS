function [trackInfo,errFlag] = evalGapsMergesSplits(tracksFinal)
%EVALGAPSMERGESSPLITS calculates displacements over gaps, merges and splits and compares them to linking displacements
%
%SYNOPSIS
%
%INPUT
%
%OUTPUT
%
%Khuloud Jaqaman, March 2010

%% Output

trackInfo = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--evalGapsMergesSplits: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%% Gap evaluation

%convert tracks from structure format to matrix format
tracksMat = convStruct2MatIgnoreMS(tracksFinal);
numTracksMat = size(tracksMat,1);

%initialize output variable
trackInfo = repmat(struct('dispGap',[]),numTracksMat,1);

%find track gaps
gapInfo = findTrackGaps(tracksMat);

%get track start and end times
trackSEL = getTrackSEL(tracksMat);

%go over all tracks
for iTrack = 1 : numTracksMat
    
    %get track coordinates
    coordTrack(:,1) = tracksMat(iTrack,1:8:end)';
    coordTrack(:,2) = tracksMat(iTrack,2:8:end)';
    coordTrack(:,3) = tracksMat(iTrack,3:8:end)';
    
    %find gaps belonging to this track
    gapsCurrentTrack = gapInfo(gapInfo(:,1)==iTrack,:);
    numGaps = size(gapsCurrentTrack,1);
    
    %initializen variable for output
    dispGap = NaN(numGaps,4);
    
    if numGaps ~= 0
        
        %get track start times and lengths
        gapStartTime = gapsCurrentTrack(:,3);
        gapLength = gapsCurrentTrack(:,4);
        
        %get unique gap lengths
        gapLengthUnique = unique(gapLength);
        numUniqueGapLengths = length(gapLengthUnique);
        
        %chop track into parts that don't contain gaps
        %do this by replicating track coordinates numGap+1 times and
        %retaining in each replica only the coordinates of one non-gap part
        %e.g. if there are 2 gaps then there will be 3 parts containing the
        %coordinates from beginning till 1st gap, between 1st gap and 2nd gap,
        %and from after 2nd gap till end
        coordTrackChopped = repmat(coordTrack,[1 1 numGaps+1]);
        coordTrackChopped(gapStartTime(1):end,:,1) = NaN;
        for iGap = 2 : numGaps
            coordTrackChopped(1:gapStartTime(iGap-1),:,iGap) = NaN;
            coordTrackChopped(gapStartTime(iGap):end,:,iGap) = NaN;
        end
        coordTrackChopped(1:gapStartTime(end),:,end) = NaN;
        
        %find the rows where each non-gap track part starts and ends
        rowStart = [trackSEL(iTrack,1); gapStartTime+gapLength];
        rowEnd = [gapStartTime-1; trackSEL(iTrack,2)];
        
        %calculate displacements over gaps
        for iGap = 1 : numGaps
            dispGap(iGap,2) = sqrt( sum( ( ...
                coordTrackChopped(rowStart(iGap+1),:,iGap+1) - ...
                coordTrackChopped(rowEnd(iGap),:,iGap) ).^2 ) );
        end
        
        %put gap length as first column of dispGap
        dispGap(:,1) = gapLength;
        
        %calculate average displacements in non-gap parts over same time
        %intervals
        for iGapLength = 1 : numUniqueGapLengths
            
            %get gap length of interest
            gapLengthCurrent = gapLengthUnique(iGapLength);
            
            %actual interval of disappearance is gapLength + 1
            disappearInt = gapLengthCurrent + 1;
            
            %get displacements over this disappearance interval
            dispNonGapCurrent = sqrt( sum( ( coordTrackChopped(disappearInt+1:end,:,:) - ...
                coordTrackChopped(1:end-disappearInt,:,:) ).^2 ,2) );
            dispNonGapCurrent = dispNonGapCurrent(:);
            dispNonGapCurrent = dispNonGapCurrent(~isnan(dispNonGapCurrent));
            numObs = length(dispNonGapCurrent);
            
            %calculate average displacement over this time interval
            aveDispNonGapCurrent = mean(dispNonGapCurrent);
            
            %store average displacement and number of observations it is
            %calculated from to compare to gap displacements
            indxThisGapLength = find(gapLength==gapLengthCurrent);
            dispGap(indxThisGapLength,3) = aveDispNonGapCurrent;
            dispGap(indxThisGapLength,4) = numObs;
            
        end
        
    end
    
    %store this track's gap displacement information in structure arrawy
    %for output
    trackInfo(iTrack).dispGap = dispGap;
    
end
