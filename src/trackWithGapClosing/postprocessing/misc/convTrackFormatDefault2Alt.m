function compTracksAlt = convTrackFormatDefault2Alt(compTracks)
%CONVTRACKFORMATDEFAULT2ALT converts compound tracks from default format to alternative format
%
%SYNOPSIS compTracksAlt = convTrackFormatDefault2Alt(compTracks)
%
%INPUT  compTracks   : Compound tracks, in the format of tracksFinal as
%                      output by trackCloseGapsKalman.
%
%OUTPUT compTracksAlt: Compound tracks in alternative format, where a
%                      merge/split leads to new track segments, without
%                      continuation of the merging/splitting segments
%                      
%REMARKS Code does not handle properly the case of multiple merges with or
%        multiple splits from the same track segment at the same time. Bug
%        must be fixed. This situation arises with simulated data.
%        Idea: Bug can be fixed by grouping events together if they are
%        simultaneous, instead of going over them sequentially.
%
%Khuloud Jaqaman, February 2009

%% Output

%initialize output
compTracksAlt = compTracks;

%% Input

%get number of compound tracks
numTracks = length(compTracks);

%% Conversion

%go over all compound tracks
for iTrack = 1 : numTracks

    %get this compound track's information
    seqOfEvents = compTracks(iTrack).seqOfEvents;
    tracksFeatIndx = compTracks(iTrack).tracksFeatIndxCG;
    tracksCoordAmp = compTracks(iTrack).tracksCoordAmpCG;

    %get unique event times
    eventTimes = unique(seqOfEvents(:,1));

    %go over unique event times
    doubleFreq = 0;
    for iEvent = 1 : length(eventTimes)

        %find events happening at this event time
        indxEventsAtEventTime = find(seqOfEvents(:,1)==eventTimes(iEvent));
        numEventsAtEventTime = length(indxEventsAtEventTime);

        %if there is more than 1 (possibly indicating both a merge and a
        %split)
        if numEventsAtEventTime > 1

            %sort events such that splits come before merges
            seqOfEventsTmp = sortrows(seqOfEvents(indxEventsAtEventTime,:),2);
            seqOfEvents(indxEventsAtEventTime,:) = seqOfEventsTmp;

            %if there are segments which exhibit a split and a merge in the
            %same frame, raise flag that this compound track will have
            %doubled sampling frequency
            probSegments = intersect(...
                seqOfEventsTmp(seqOfEventsTmp(:,2)==1&~isnan(seqOfEventsTmp(:,4)),4),...
                seqOfEventsTmp(seqOfEventsTmp(:,2)==2&~isnan(seqOfEventsTmp(:,4)),4));
            if ~isempty(probSegments)
                doubleFreq = 1;
            end
            
        end %(if numEventsAtEventTime > 1)

    end %(for iEvent = 1 : length(eventTimes))
    
    %copy information for alternative format, with sampling frequency
    %doubling if necessary
    if doubleFreq
        seqOfEventsAlt = seqOfEvents;
        seqOfEventsAlt(seqOfEventsAlt(:,2)==1,1) = ...
            seqOfEventsAlt(seqOfEventsAlt(:,2)==1,1) - 0.5;
        [numRows,numCols] = size(tracksFeatIndx);
        tracksFeatIndxAlt = zeros(numRows,2*numCols);
        tracksFeatIndxAlt(:,1:2:end) = tracksFeatIndx;
        tracksFeatIndxAlt(:,2:2:end) = tracksFeatIndx;
        tracksCoordAmpAlt = NaN(numRows,16*numCols);
        
        %090214 - Robel Yirdaw
        %Convert sparse matrix back to full form to avoid bottleneck at
        %lines 169 - 172.        
        if (issparse(tracksCoordAmp))
            tracksCoordAmp = full(tracksCoordAmp);
        end
        tracksCoordAmpAlt(:,1:16:end)  = tracksCoordAmp(:,1:8:end);
        tracksCoordAmpAlt(:,2:16:end)  = tracksCoordAmp(:,2:8:end);
        tracksCoordAmpAlt(:,3:16:end)  = tracksCoordAmp(:,3:8:end);
        tracksCoordAmpAlt(:,4:16:end)  = tracksCoordAmp(:,4:8:end);
        tracksCoordAmpAlt(:,5:16:end)  = tracksCoordAmp(:,5:8:end);
        tracksCoordAmpAlt(:,6:16:end)  = tracksCoordAmp(:,6:8:end);
        tracksCoordAmpAlt(:,7:16:end)  = tracksCoordAmp(:,7:8:end);
        tracksCoordAmpAlt(:,8:16:end)  = tracksCoordAmp(:,8:8:end);
        tracksCoordAmpAlt(:,9:16:end)  = tracksCoordAmp(:,1:8:end);
        tracksCoordAmpAlt(:,10:16:end) = tracksCoordAmp(:,2:8:end);
        tracksCoordAmpAlt(:,11:16:end) = tracksCoordAmp(:,3:8:end);
        tracksCoordAmpAlt(:,12:16:end) = tracksCoordAmp(:,4:8:end);
        tracksCoordAmpAlt(:,13:16:end) = tracksCoordAmp(:,5:8:end);
        tracksCoordAmpAlt(:,14:16:end) = tracksCoordAmp(:,6:8:end);
        tracksCoordAmpAlt(:,15:16:end) = tracksCoordAmp(:,7:8:end);
        tracksCoordAmpAlt(:,16:16:end) = tracksCoordAmp(:,8:8:end);
    else
        seqOfEventsAlt = seqOfEvents;
        tracksFeatIndxAlt = tracksFeatIndx;
        %082114 - Robel Yirdaw
        %Convert sparse matrix back to full form to avoid bottleneck at
        %lines 169 - 172.
        if (issparse(tracksCoordAmp))
            tracksCoordAmpAlt = full(tracksCoordAmp);
        else
            tracksCoordAmpAlt = tracksCoordAmp;
        end
    end
    
    %shift time in seqOfEventsAlt to make the track start at frame 1 or 0.5
    %if doubling sampling frequency
    seqOfEventsAlt(:,1) = seqOfEventsAlt(:,1) - seqOfEvents(1,1) + 1;
    
    %get number of segments and events in default format
    numSegmentsDef = size(tracksFeatIndx,1);
    numSegmentsAlt = numSegmentsDef;
    numEventsDef = size(seqOfEvents,1);
    
    %find all merging and splitting events
    msEvents = find(~isnan(seqOfEvents(:,4)));
    numMSEvents = length(msEvents);

    %reserve memory for matrix storing segment correspondence
    alt2defSegCorr = zeros(numMSEvents,2);
    
    %transpose matrices to speed up calculations
    tracksFeatIndxAlt = tracksFeatIndxAlt';
    tracksCoordAmpAlt = tracksCoordAmpAlt';
    
    %go over merging and splitting events
    for iEventTmp = 1 : numMSEvents
        %get event index
        iEvent = msEvents(iEventTmp);
        
        %get time and type of event, and "original" segment that got merged
        %with or split from
        eventTime = seqOfEventsAlt(iEvent,1);
        eventType = seqOfEventsAlt(iEvent,2);
        originalSegment = seqOfEventsAlt(iEvent,4);

        %add 1 to number of segments
        numSegmentsAlt = numSegmentsAlt + 1;
        
        %store segment correspondence
        segmentCorrespond = originalSegment;
        while segmentCorrespond > numSegmentsDef
            segmentCorrespond = alt2defSegCorr(alt2defSegCorr(:,1)==segmentCorrespond,2);
        end
        alt2defSegCorr(iEventTmp,:) = [numSegmentsAlt segmentCorrespond];
        
        %separate chunk after merging or splitting from original segment ...
        
        %update feature connectivity matrix
        %         tracksFeatIndxAlt(numSegmentsAlt,(doubleFreq+1)*eventTime:end) = tracksFeatIndxAlt(...
        %             originalSegment,(doubleFreq+1)*eventTime:end);
        %         tracksFeatIndxAlt(numSegmentsAlt,1:(doubleFreq+1)*eventTime-1) = 0;
        %         tracksFeatIndxAlt(originalSegment,(doubleFreq+1)*eventTime:end) = 0;
        tracksFeatIndxAlt((doubleFreq+1)*eventTime:end,numSegmentsAlt) = tracksFeatIndxAlt(...
            (doubleFreq+1)*eventTime:end,originalSegment);
        tracksFeatIndxAlt(1:(doubleFreq+1)*eventTime-1,numSegmentsAlt) = 0;
        tracksFeatIndxAlt((doubleFreq+1)*eventTime:end,originalSegment) = 0;
        
        %update matrix of coordinates and amplitudes
        %         tracksCoordAmpAlt(numSegmentsAlt,8*((doubleFreq+1)*eventTime-1)+1:end) = ...
        %             tracksCoordAmpAlt(originalSegment,8*((doubleFreq+1)*eventTime-1)+1:end);
        %         tracksCoordAmpAlt(numSegmentsAlt,1:8*((doubleFreq+1)*eventTime-1)) = NaN;
        %         tracksCoordAmpAlt(originalSegment,8*((doubleFreq+1)*eventTime-1)+1:end) = NaN;
        tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)+1:end,numSegmentsAlt) = ...
            tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)+1:end,originalSegment);
        tracksCoordAmpAlt(1:8*((doubleFreq+1)*eventTime-1),numSegmentsAlt) = NaN;
        tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)+1:end,originalSegment) = NaN;
        
        %update sequence of events based on event type
        switch eventType

            case 1 %if it's a split
                
                %add new event for splitting of original segment into the
                %new separated segment
                seqOfEventsAlt(end+1,:) = [eventTime eventType numSegmentsAlt ...
                    originalSegment];

            case 2 %if it's a merge

                %update index of segment merged into
                seqOfEventsAlt(iEvent,4) = numSegmentsAlt;
                
                %add new event for merging of original segment into the new
                %separated segment
                seqOfEventsAlt(end+1,:) = [eventTime eventType ...
                    originalSegment numSegmentsAlt];
                
                %update feature connectivity and coordinate matrices of
                %sampling frequency is doubled
                if doubleFreq

                    %get index of merging segment
                    mergingSegment = seqOfEventsAlt(iEvent,3);

                    %update its feature connectivity matrix and matrix of
                    %coordinates and amplitudes
                    %                     tracksFeatIndxAlt(mergingSegment,(doubleFreq+1)*eventTime-1) = ...
                    %                         tracksFeatIndxAlt(mergingSegment,(doubleFreq+1)*eventTime-2);
                    %                     tracksCoordAmpAlt(mergingSegment,8*((doubleFreq+1)*eventTime-1)-7:...
                    %                         8*((doubleFreq+1)*eventTime-1)) = tracksCoordAmpAlt(...
                    %                         mergingSegment,8*((doubleFreq+1)*eventTime-1)-15:...
                    %                         8*((doubleFreq+1)*eventTime-1)-8);
                    tracksFeatIndxAlt((doubleFreq+1)*eventTime-1,mergingSegment) = ...
                        tracksFeatIndxAlt((doubleFreq+1)*eventTime-2,mergingSegment);
                    tracksCoordAmpAlt(8*((doubleFreq+1)*eventTime-1)-7:...
                        8*((doubleFreq+1)*eventTime-1),mergingSegment) = tracksCoordAmpAlt(...
                        8*((doubleFreq+1)*eventTime-1)-15:...
                        8*((doubleFreq+1)*eventTime-1)-8,mergingSegment);

                end

        end

        %if the original segment is involved in later events,
        %replace it with the new additional segment
        indx = find(seqOfEventsAlt(iEvent+1:numEventsDef,3)==originalSegment);
        seqOfEventsAlt(iEvent+indx,3) = numSegmentsAlt;
        indx = find(seqOfEventsAlt(iEvent+1:numEventsDef,4)==originalSegment);
        seqOfEventsAlt(iEvent+indx,4) = numSegmentsAlt;

    end

    %transpose matrices back
    tracksFeatIndxAlt = tracksFeatIndxAlt';
    tracksCoordAmpAlt = tracksCoordAmpAlt';
    
    %sort sequence of events in ascending chronological order
    seqOfEventsAlt = sortrows(seqOfEventsAlt,1);
    
    %get unique event times
    eventTimes = unique(seqOfEventsAlt(:,1));
    
    %go over unique event times
    for iEvent = 1 : length(eventTimes)
        
        %find events happening at this event time
        indxEventsAtEventTime = find(seqOfEventsAlt(:,1)==eventTimes(iEvent));
        numEventsAtEventTime = length(indxEventsAtEventTime);

        %if there are more than 2 [1=initiation or termination (OK case);
        %2 = merge or split (OK case); 3 = merge or split AND initiation or
        %termination (possible "problem" case), etc. (all possible
        %"problem" cases)]
        if numEventsAtEventTime > 2
            
            %sort events such that splits/starts come before merges/ends
            seqOfEventsTmp = sortrows(seqOfEventsAlt(indxEventsAtEventTime,:),2);
            
            %extract rows documenting starts, splits, ends and merges
            seqOfEventsTmpStart = seqOfEventsTmp(seqOfEventsTmp(:,2)==1&...
                isnan(seqOfEventsTmp(:,4)),:);
            seqOfEventsTmpSplit = seqOfEventsTmp(seqOfEventsTmp(:,2)==1&...
                ~isnan(seqOfEventsTmp(:,4)),:);
            seqOfEventsTmpEnd = seqOfEventsTmp(seqOfEventsTmp(:,2)==2&...
                isnan(seqOfEventsTmp(:,4)),:);
            seqOfEventsTmpMerge = seqOfEventsTmp(seqOfEventsTmp(:,2)==2&...
                ~isnan(seqOfEventsTmp(:,4)),:);

            %make sure that every merge/split is documented on two
            %consecutive rows
            seqOfEventsTmpSplit = sortrows(seqOfEventsTmpSplit,4);
            seqOfEventsTmpMerge = sortrows(seqOfEventsTmpMerge,4);

            %put documentation back in sequence of events
            seqOfEventsAlt(indxEventsAtEventTime,:) = [seqOfEventsTmpStart; ...
                seqOfEventsTmpSplit; seqOfEventsTmpMerge; seqOfEventsTmpEnd];
            
        end
        
    end
    
    %add back time offset to get real times
    seqOfEventsAlt(:,1) = seqOfEventsAlt(:,1) + seqOfEvents(1,1) - 1;
    
    %store in output structure
    %082014 - Robel Yirdaw
    %Save tracksCoordAmpCG as a sparse matrix if originally sparse    
    compTracksAlt(iTrack).seqOfEvents = seqOfEventsAlt;
    compTracksAlt(iTrack).tracksFeatIndxCG = tracksFeatIndxAlt;
    if (issparse(tracksCoordAmp))
        tracksCoordAmpAlt(isnan(tracksCoordAmpAlt)) = 0;
        compTracksAlt(iTrack).tracksCoordAmpCG = sparse(tracksCoordAmpAlt);
    else
        compTracksAlt(iTrack).tracksCoordAmpCG = tracksCoordAmpAlt;
    end
    compTracksAlt(iTrack).alt2defSegmentCorrespond = alt2defSegCorr;
    
end

%% ~~~ the end ~~~
