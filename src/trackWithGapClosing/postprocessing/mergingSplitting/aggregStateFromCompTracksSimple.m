function compTracks = aggregStateFromCompTracksSimple(compTracks,intensityInfo)
%AGGREGSTATEFROMCOMPTRACKS recovers particle aggregation states from compound tracks
%
%SYNOPSIS compTracks = aggregStateFromCompTracks(compTracks,intensityInfo)
%
%INPUT  compTracks   : Compound tracks, in the format of tracksFinal as
%                      output by trackCloseGapsKalman.
%       intensityInfo: Row vector with unit intensity mean and standard
%                      deviation (e.g. the intensity of a single
%                      fluorophore labeling a single receptor).
%
%OUTPUT compTracks   : Same as input but with additional field
%           .aggregState: Intenger number of "units" (e.g. receptors)
%                         within each detected particle (or spot).

%Khuloud Jaqaman, February 2009

%% Input

%assign default intensityInfo if not input
if nargin < 2 || isempty(intensityInfo)
    intensityInfo = [];
end

%get number of compound tracks
numTracks = length(compTracks);

%% Calculation

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get this compound track's information
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    tracksFeatIndx = compTracks(iTrack).tracksFeatIndxCG;
    tracksAmp = compTracks(iTrack).tracksCoordAmpCG(:,4:8:end);

    %shift time in seqOfEvents to make the track start at frame 1
    seqOfEvents(:,1) = seqOfEvents(:,1) - seqOfEvents(1,1) + 1;
    
    %initialize matrix storing aggregation state
    if isempty(intensityInfo) %if no information is given on the unit intensity
        aggregStateMat = tracksFeatIndx; %initialize every branch with 1 unit
        aggregStateMat = double((aggregStateMat ~= 0));
        aggregStateMat(aggregStateMat==0) = NaN;
    else %if there is unit intensity information
        
    end
    
    %find all merging and splitting events
    msEvents = find(~isnan(seqOfEvents(:,4)));
    
    %go over merging and splitting events and modify aggregation state
    %accordingly
    for iEvent = msEvents'
        
        %get nature and time of event
        eventTime = seqOfEvents(iEvent,1);
        eventType = seqOfEvents(iEvent,2);
        
        %update aggregation state based on event type
        switch eventType
            
            case 1 %split
                
                %find the original segment that the new segment split from
                segmentOriginal = seqOfEvents(iEvent,4);
                
                %if the original segment has only 1 molecule ...
                if aggregStateMat(segmentOriginal,eventTime-1) == 1

                    %add 1 molecule to its aggregation state before the
                    %split
                    aggregStateMat(segmentOriginal,1:eventTime-1) = ...
                        aggregStateMat(segmentOriginal,1:eventTime-1) + 1;

                    %determine whether the original segment itself splitted
                    %from another segment
                    origSegmentStartEvent = find(seqOfEvents(:,2)==1 & ...
                        seqOfEvents(:,3)==segmentOriginal);
                    splitting = ~isnan(seqOfEvents(origSegmentStartEvent,4));

                    %as long as an "original" segment splitted from an "original
                    %original" segment, loop and update aggregation state
                    %stop when one segment finally starts by an appearance
                    while splitting

                        %find split time and segment it split from
                        splitTime = seqOfEvents(origSegmentStartEvent,1);
                        segmentOrigOrig = seqOfEvents(origSegmentStartEvent,4);

                        %update the aggregation state of the segment it split
                        %from, before the split
                        aggregStateMat(segmentOrigOrig,1:splitTime-1) = ...
                            aggregStateMat(segmentOrigOrig,1:splitTime-1) + 1;

                        %determine whether this "original original" segment
                        %splitted from another segment
                        origSegmentStartEvent = find(seqOfEvents(:,2)==1 & ...
                            seqOfEvents(:,3)==segmentOrigOrig);
                        splitting = ~isnan(seqOfEvents(origSegmentStartEvent,4));

                    end
                    
                else %if the original segment has more than 1 molecule ...
                    
                    %subtract 1 molecule from its aggregation state after
                    %the split
                    aggregStateMat(segmentOriginal,eventTime:end) = ...
                        aggregStateMat(segmentOriginal,eventTime:end) - 1;
                    
                end

            case 2 %merge
                
                %find the two merging segments
                segmentOriginal = seqOfEvents(iEvent,4);
                segmentMerging = seqOfEvents(iEvent,3);
                
                %update its aggregation state as the sum of the aggregation
                %states before merging
                aggregStateMat(segmentOriginal,eventTime:end) = ...
                    aggregStateMat(segmentOriginal,eventTime-1) + ...
                    aggregStateMat(segmentMerging,eventTime-1) + ...
                    aggregStateMat(segmentOriginal,eventTime:end) - ...
                    aggregStateMat(segmentOriginal,eventTime:end);
                
        end %(switch eventType)
        
    end %(for iEvent = msEvents')
    
    compTracks(iTrack).aggregState = aggregStateMat;
    
end %(for iTrack = 1 : numTracks)

%% ~~~ the end ~~~