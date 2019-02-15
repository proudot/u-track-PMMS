function [msStats,msTimeInfo,mergesInfo,splitsInfo] = calcStatsMSSimple(...
    tracks,minTrackLen,probDim)
%MSSTATS calculate some merge/split statistics
%
%SYNOPSIS [msStats,msTimeInfo,mergesInfo,splitsInfo] = calcStatsMSSimple(...
%    tracks,minTrackLen,probDim)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%
%OUTPUT msStats    : Row vector with entries: 
%                    1st: number of features per frame;.
%                    2nd/3rd: Number of merges/splits per feature.
%                    4th: Total number of tracks (length >= minTrackLen).
%                    5th: Fraction of tracks that are linear.
%                    6th: Fraction of features undergoing linear motion.
%                    7th/8th: Fraction of merges/splits in linear tracks.
%                    9th/10th: Probability of a feature merging/splitting
%                              while in a linear track.
%                    11th/12th: Probability of a feature merging/splitting
%                               while in a Brownian track.
%                    13th/14th: Ratio of linear probability to Brownian
%                               probability for merges/splits.
%       msTimeInfo : Structure with field 'brown' and 'linear' for Brownian
%                    and linear tracks. Each field is a structure
%                    containing the fields:
%           .numTracks      : 4 entries: # of tracks with merges and
%                             splits, # of tracks with only merges, # of
%                             tracks with only splits, and # of tracks
%                             without any merges or splits.
%           .timeMerge2Split: 2-column vector where first row is time
%                             between a merge and a consecutive split and
%                             second column is the weight of that
%                             observation.
%           .timeSplit2Merge: 2-column vector where first row is time
%                             between a split and a consecutive merge and
%                             second column is the weight of that
%                             observation.
%       mergesInfo : Output of findMergesSplits.
%       splitsInfo : Output of findMergesSplits.
%
%Khuloud Jaqaman, December 2007

%% input

if nargin < 1 || isempty(tracks)
    disp('calcStatsMSSimple: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%% features

%get average number of features per frame
numFeatTot = 0; 
for iTrack = 1 : numTracks
    numFeatTot = numFeatTot + length(find(tracks(iTrack).tracksFeatIndxCG)); 
end
aveFeatPerFrame = numFeatTot / numFrames;

%% track types

%find number of tracks per type
criteria.trackType = 0;
indxBrown = chooseTracks(tracks,criteria);
numTracksBrown = length(indxBrown);
criteria.trackType = 1;
indxLin = chooseTracks(tracks,criteria);
numTracksLin = length(indxLin);
clear criteria

%calculate fraction of tracks that are linear
fracTrackLin = numTracksLin / (numTracksBrown + numTracksLin);

%calculate number of features undergoing linear motion
numFeatLin = 0;
for iTrack = indxLin'
    numFeatLin = numFeatLin + length(find(tracks(iTrack).tracksFeatIndxCG));
end

%get fraction of features undergoing linear motion
fracFeatLin = numFeatLin / numFeatTot;

%% merges/splits vs. features

%locate merges and splits in tracks
[mergesInfo,splitsInfo] = findMergesSplits(tracks,probDim);

%get total number of merges and splits
numMergesTot = sum(mergesInfo(:,3));
numSplitsTot = sum(splitsInfo(:,3));

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate average number of merges/splits per feature
aveMergePerFeat = aveMergePerFrame / aveFeatPerFrame;
aveSplitPerFeat = aveSplitPerFrame / aveFeatPerFrame;

%% merges/splits vs. track type

%get number of merges/splits happening in linear tracks
numMergesLin = sum(mergesInfo(mergesInfo(:,2)==1,3));
numSplitsLin = sum(splitsInfo(splitsInfo(:,2)==1,3));

%get their fraction of total number
fracMergesLin = numMergesLin / numMergesTot;
fracSplitsLin = numSplitsLin / numSplitsTot;

%calculate the probability of a feature undergoing a merge/split while in a
%linear track
probMergeLin = fracMergesLin * aveMergePerFeat / fracFeatLin;
probSplitLin = fracSplitsLin * aveSplitPerFeat / fracFeatLin;

%calculate the probability of a feature undergoing a merge/split while in a
%Brownian track
probMergeBr = (1-fracMergesLin) * aveMergePerFeat / (1-fracFeatLin);
probSplitBr = (1-fracSplitsLin) * aveSplitPerFeat / (1-fracFeatLin);

%% time from merges to splits and vice versa

for iType = 0 : 1

    if iType == 1
        trackType = 'Lin';
    else
        trackType =  'Brown';
    end
    
    %initialize some variables
    numTrackNoMS = 0;
    numTrackOnlyM = 0;
    numTrackOnlyS = 0;
    numTrackMS = 0;
    timeMerge2Split = [];
    timeSplit2Merge = [];
    eval(['indxTracks = indx' trackType ';']);

    %go over all tracks of this type ...
    for iTrack = indxTracks'

        %get track's sequence of events
        seqOfEvents = tracks(iTrack).seqOfEvents;

        %find merge and split times in track
        msTime = seqOfEvents(~isnan(seqOfEvents(:,4)),[1 2]);
        msTime(msTime(:,2)==1,1) = -msTime(msTime(:,2)==1,1);
        msTime = msTime(:,1);

        %take action based on whether there are merges and splits
        if isempty(msTime) %if there are no merges and no splits

            %add one to counter of tracks without merges and splits
            numTrackNoMS = numTrackNoMS + 1;

        elseif all(msTime > 0) %if there are merges but no splits

            %add one to counter of tracks with only merges
            numTrackOnlyM = numTrackOnlyM + 1;

        elseif all(msTime < 0) %if there are splits but no merges

            %add one to counter of tracks with only splits
            numTrackOnlyS = numTrackOnlyS + 1;

        else %if there are both merges and splits

            %add one to counter of tracks with both merges and
            %splits
            numTrackMS = numTrackMS + 1;

            %get number of merges and splits to consider
            numMS = length(msTime);

            indxMS = 1;
            while indxMS < numMS

                %get index of initial event and its type
                typeEvent1 = sign(msTime(indxMS));
                indxEvent1 = indxMS;

                %keep increasing indxMS until you find an event different
                %from initial event
                typeEvent2 = typeEvent1;
                while typeEvent2 == typeEvent1 && indxMS < numMS
                    indxMS = indxMS + 1;
                    typeEvent2 = sign(msTime(indxMS));
                end
                indxEvent2 = indxMS;

                %keep increasing indxMS until you find an event again
                %similar to initial event
                typeEvent3 = typeEvent2;
                while typeEvent3 ~= typeEvent1 && indxMS < numMS
                    indxMS = indxMS + 1;
                    typeEvent3 = sign(msTime(indxMS));
                end
                if typeEvent3 ~= typeEvent2
                    indxMS = indxMS - 1;
                end
                indxEvent3 = indxMS;

                if typeEvent2 ~= typeEvent1

                    %calculate total number of combinations for calculating
                    %time between merges and splits (or vice versa)
                    numCombination = (indxEvent2-indxEvent1)*(indxEvent3-indxEvent2+1);

                    %go over all combinations and calculate time
                    indxComb = 0;
                    timeBetweenMS = zeros(numCombination,1);
                    for indx1 = 1 : indxEvent2-indxEvent1
                        for indx2 = 1 : indxEvent3-indxEvent2+1
                            indxComb = indxComb + 1;
                            timeBetweenMS(indxComb) = abs(msTime(indx2+indxEvent2-1)) - ...
                                abs(msTime(indx1+indxEvent1-1));
                        end
                    end

                    %store the times and their weights based on whether we
                    %looked at a merge to split or a split to merge
                    if typeEvent1 > 0
                        timeMerge2Split = [timeMerge2Split; [timeBetweenMS ...
                            (1/numCombination)*ones(numCombination,1)]];
                    else
                        timeSplit2Merge = [timeSplit2Merge; [timeBetweenMS ...
                            (1/numCombination)*ones(numCombination,1)]];
                    end

                end

                %update indxMS to look at next merge-to-split or
                %split-to-merge event
                indxMS = indxEvent2;

            end %(while indxMS <= numMS - 1)

        end %(if isempty(msTime) ... elseif ...)

    end %(for iTrack = 1 : indxTracks')
    
    %store track numbers and distributions
    eval(['numTracks' trackType ' = [numTrackMS numTrackOnlyM numTrackOnlyS numTrackNoMS];'])
    eval(['timeMerge2Split' trackType ' = timeMerge2Split;'])
    eval(['timeSplit2Merge' trackType ' = timeSplit2Merge;'])
    
end
    
%% output

msStats = [aveFeatPerFrame aveMergePerFeat aveSplitPerFeat ...
    numTracks fracTrackLin fracFeatLin fracMergesLin fracSplitsLin ...
    probMergeLin probSplitLin probMergeBr probSplitBr ...
    probMergeLin./probMergeBr probSplitLin./probSplitBr];

msTimeInfo.brown.numTracks = numTracksBrown;
msTimeInfo.brown.timeMerge2Split = timeMerge2SplitBrown;
msTimeInfo.brown.timeSplit2Merge = timeSplit2MergeBrown;
msTimeInfo.linear.numTracks = numTracksLin;
msTimeInfo.linear.timeMerge2Split = timeMerge2SplitLin;
msTimeInfo.linear.timeSplit2Merge = timeSplit2MergeLin;

%% ~~~ the end ~~~



