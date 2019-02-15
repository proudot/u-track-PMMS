function [segmentSetOUT,segmentErr,callCount] = processSegmentEvents(prevEventTime,segID,segmentSetIN,...
    segmentErr,seqOfEvents,tracksAmp,intensityInfo,callCount)
%PROCESSSEGMENTEVENTS recursively goes through sequence of events and 
%processes merging and splitting events constructing a detailed 
%relationship of the segments.
%
%   INPUT:
%           prevEventTime:      the previous iteration point
%           segID:              ID of segment currently being processed
%           segmentSetIN:       table with each row corresponding to a
%                               segment and the following columns:
%                               col. 1: ID of parent segment 1
%                               col. 2: ID of parent segment 2, if any
%                               col. 3: ID of child segment 1, if any
%                               col. 4: ID of child segment 2, if any
%                               col. 5: ID of sibling segment, if any
%                               col. 6: ID of partner segment, if any
%                               col. 7: intensity amplitude of segID
%                               col. 8: segID's size from intensity
%                               col. 9: segID's size from events
%                               col. 10: a boolean indicating whether a
%                                        size repair is being performed
%                                        (used by the repair functions).
%           segmentErr:         a 2 element cell vector containing segment
%                               IDs with intensity size errors in element 1
%                               and IDs with event size errors in element 2
%           seqOfEvents:        sequence of events from compTracks in
%                               alternative format
%           tracksAmp:          intensity amplitudes from compTracks
%           intensityInfo:      intensity quantum used in the simulation,
%                               in the same format, i.e. [1 0.3]
%           callCount:          tracks number of recursive calls
%
%   OUTPUT:
%           segementSetOUT:     completed table as detailed above
%           segmentErr:         collection of segment ID with errors in
%                               sizes from intensity and evetns, as
%                               described above
%           callCount:          number of recursive calls
%
%   
%   Robel Yirdaw, November 2013
%       Modified, 08/20/14
%



    segmentSetOUT = segmentSetIN;
        
    %if ((isnan(segmentSetIN(segID,1))) || isnan(segmentSetIN(segID,3)) )  
        %segID is the parent segment of an event (split or merge). In the
        %case of a split it will be on the 4th column of sequence of events
        %and in the case of a merge it will be on the 3rd column. Determine
        %the row for the next event of segID accordingly.
        if (isnan(prevEventTime))
            nextEventRow = find(((seqOfEvents(:,3) == segID) | (seqOfEvents(:,4) == segID)) & ...
                ~isnan(seqOfEvents(:,4)));
        else        
            nextEventRow = find((seqOfEvents(:,1) > prevEventTime) & ...
                ((seqOfEvents(:,3) == segID) | (seqOfEvents(:,4) == segID)) & ...
                ~isnan(seqOfEvents(:,4)));
        end

        if (~isempty(nextEventRow))
            nextEventTime = seqOfEvents(nextEventRow(1),1);
            nextEventType = seqOfEvents(nextEventRow(1),2);

            if (nextEventType == 1)
                %SPLIT EVENT
                parentID = segID;            
                child1ID = seqOfEvents(nextEventRow(1),3);
                child2ID = seqOfEvents(nextEventRow(2),3);
                %Update parent
                segmentSetOUT(parentID,3) = child1ID;
                segmentSetOUT(parentID,4) = child2ID;
                %If parent is initial segment, intensity amplitude
                %needs to be set here.
                if (isnan(prevEventTime))                
                    segmentSetOUT(parentID,7) = ...
                        mean( tracksAmp(parentID,(~isnan(tracksAmp(parentID,:)) & tracksAmp(parentID,:) ~= 0)) )/intensityInfo(1);
                    
                    segmentSetOUT(parentID,8) = round(segmentSetOUT(parentID,7));
                end                      
                %Update children
                segmentSetOUT(child1ID,1) = parentID;
                segmentSetOUT(child1ID,5) = child2ID;
                %Update intensity & size only if not previously updated.
                %Note that a segment is assigned child information after
                %all other fields have been set, when the segment is a
                %parent segment. Intensity & size information could have
                %been set when another related event was processed and so
                %it shouldn't be reset here.
                if (isnan(segmentSetOUT(child1ID,7)))
                    segmentSetOUT(child1ID,7) = ...
                        mean( tracksAmp(child1ID,(~isnan(tracksAmp(child1ID,:)) & tracksAmp(child1ID,:) ~= 0)) )/intensityInfo(1);                    
                    segmentSetOUT(child1ID,8) = round(segmentSetOUT(child1ID,7));
                end
                segmentSetOUT(child2ID,1) = parentID;
                segmentSetOUT(child2ID,5) = child1ID;
                if (isnan(segmentSetOUT(child2ID,7)))
                    segmentSetOUT(child2ID,7) = ...
                        mean( tracksAmp(child2ID,(~isnan(tracksAmp(child2ID,:)) & tracksAmp(child2ID,:) ~= 0)) )/intensityInfo(1);
                    segmentSetOUT(child2ID,8) = round(segmentSetOUT(child2ID,7));
                end

                %Check intensity amplitude values
                if ( ( (segmentSetOUT(parentID,8) == 0) || (segmentSetOUT(child1ID,8) == 0) || ...
                        (segmentSetOUT(child2ID,8) == 0) ) && (~any(segmentErr{1} == parentID)) )
                    %There is an intensity size discrepancy. Flag parent.
                    segmentErr{1}(sum(~isnan(segmentErr{1}))+1) = parentID;                    
                    %If parent has size 0, it must be set to size 1.
                    if (segmentSetOUT(parentID,8) == 0)
                        segmentSetOUT(parentID,8) = 1;
                    end
                    if ((segmentSetOUT(child1ID,8) == 0) && isnan(segmentSetOUT(child1ID,3)) )
                        segmentSetOUT(child1ID,8) = 1;
                    end
                    if ((segmentSetOUT(child2ID,8) == 0) && isnan(segmentSetOUT(child2ID,3)) )
                        segmentSetOUT(child2ID,8) = 1;
                    end
                elseif ( (segmentSetOUT(child1ID,8) > 1) && (segmentSetOUT(child2ID,8) > 1) &&...
                        (~any(segmentErr{1} == parentID)))
                    %Child segments can't both be > 1.
                    segmentErr{1}(sum(~isnan(segmentErr{1}))+1) = parentID;    
                elseif ( (segmentSetOUT(parentID,8) ~= segmentSetOUT(child1ID,8) + ...
                        segmentSetOUT(child2ID,8) ) && (~any(segmentErr{1} == parentID)) )
                    %There is an intensity size discrepancy.
                    segmentErr{1}(sum(~isnan(segmentErr{1}))+1) = parentID;
                end

                %fprintf('\n(%d) Time: %d, (%d) -> %d + %d.',callCount,nextEventTime,parentID,child1ID,child2ID);
                
                if (isnan(segmentSetIN(child1ID,3)))
                    callCount = callCount + 1;
                    [segmentSetOUT,segmentErr,callCount] = processSegmentEvents(nextEventTime,child1ID,segmentSetOUT,...
                        segmentErr,seqOfEvents,tracksAmp,intensityInfo,callCount);
                    callCount = callCount - 1;
                end
                if (isnan(segmentSetIN(child2ID,3)))
                    callCount = callCount + 1;
                    [segmentSetOUT,segmentErr,callCount] = processSegmentEvents(nextEventTime,child2ID,segmentSetOUT,...
                        segmentErr,seqOfEvents,tracksAmp,intensityInfo,callCount);     
                    callCount = callCount - 1;
                end
                %Check amplitude from events
                if ((segmentSetOUT(parentID,9) > 1) && ...
                        (segmentSetOUT(parentID,9) ~= segmentSetOUT(child1ID,9) + ...
                        segmentSetOUT(child2ID,9)) )
                    %There is an event size discrepancy.
                    %This is second condition in split block of
                    %original code (if number of molecules > 1)
                    if (~any(segmentErr{2}(:,1) == parentID))
                        currentIndx = sum(~isnan(segmentErr{2}(:,1)))+1;
                        segmentErr{2}(currentIndx,1) = parentID;
                        %111213 - not working all the time
                        segmentErr{2}(currentIndx,2) = max(segmentSetOUT(parentID,9),...
                            segmentSetOUT(child1ID,9) + segmentSetOUT(child2ID,9));
                    end
                else                                       
                    %Update parent's amplitude from events
                    segmentSetOUT(segID,9) = segmentSetOUT(child1ID,9) + ...
                        segmentSetOUT(child2ID,9);
                end

            elseif (nextEventType == 2)
                %MERGE EVENT
                %There are two parents. Incoming segment is one.
                parent1ID = segID;
                parent2Row = ((seqOfEvents(:,1) == seqOfEvents(nextEventRow,1) ) & ...
                     (seqOfEvents(:,2) == 2) & (seqOfEvents(:,3) ~= parent1ID) & ...
                     (seqOfEvents(:,4) == seqOfEvents(nextEventRow,4)) );
                parent2ID = seqOfEvents(parent2Row,3);
                childID = seqOfEvents(nextEventRow(1),4);        

                %fprintf('\n (%d) Time: %d, (%d) + %d -> %d.',callCount,nextEventTime,parent1ID,parent2ID,childID);
                %Update parents
                if (isnan(segmentSetIN(parent1ID,3)))
                    segmentSetOUT(parent1ID,3) = childID;
                    segmentSetOUT(parent1ID,6) = parent2ID;
                    %segmentSetOUT(parent2ID,3) = childID;
                    %segmentSetOUT(parent2ID,6) = parent1ID;
                    %If merge of initial segments, update their intensities.
                    if (isnan(segmentSetOUT(parent1ID,7)))
                        segmentSetOUT(parent1ID,7) =...
                            mean( tracksAmp(parent1ID,(~isnan(tracksAmp(parent1ID,:)) & tracksAmp(parent1ID,:) ~= 0)) )/intensityInfo(1);
                        segmentSetOUT(parent1ID,8) = round(segmentSetOUT(parent1ID,7));
                    end
                    
                    if (isnan(segmentSetOUT(parent2ID,7)))
                        segmentSetOUT(parent2ID,7) =...
                            mean( tracksAmp(parent2ID,(~isnan(tracksAmp(parent2ID,:)) & tracksAmp(parent2ID,:) ~= 0)) )/intensityInfo(1);                        
                        segmentSetOUT(parent2ID,8) = round(segmentSetOUT(parent2ID,7));
                    end
                        
                    %Update child
                    segmentSetOUT(childID,1) = parent1ID;
                    segmentSetOUT(childID,2) = parent2ID;
                    if  (isnan(segmentSetOUT(childID,7)))
                        segmentSetOUT(childID,7) = ...
                            mean( tracksAmp(childID,(~isnan(tracksAmp(childID,:)) & tracksAmp(childID,:) ~= 0)) )/intensityInfo(1);
                        segmentSetOUT(childID,8) = round(segmentSetOUT(childID,7));   
                    end
                    
                    %Check intensity amplitude values
                    if ( ( (segmentSetOUT(parent1ID,8) == 0) || (segmentSetOUT(parent2ID,8) == 0) ||...
                            (segmentSetOUT(childID,8) == 0) ) && ... 
                            (~any(segmentErr{1} == parent1ID) && ~any(segmentErr{1} == parent2ID)) )
                        %There is an intensity size discrepancy.
                        segmentErr{1}(sum(~isnan(segmentErr{1}))+1) = parent1ID;   
                        %If either parent has size 0, it must be set to size 1.                        
                        if (segmentSetOUT(parent1ID,8) == 0)
                            segmentSetOUT(parent1ID,8) = 1;
                        end
                        if (segmentSetOUT(parent2ID,8) == 0)
                            segmentSetOUT(parent2ID,8) = 1;
                        end
                        if ((segmentSetOUT(childID,8) == 0) && isnan(segmentSetOUT(childID,3)) )
                            segmentSetOUT(childID,8) = 1;
                        end
                    elseif ( (segmentSetOUT(parent1ID,8) > 1) && (segmentSetOUT(parent2ID,8) > 1) &&...
                             ~any(segmentErr{1} == parent1ID) && ~any(segmentErr{1} == parent2ID) )
                        %Both parents can't have size > 1.
                        segmentErr{1}(sum(~isnan(segmentErr{1}))+1) = parent1ID;                        
                        
                    elseif ( ( segmentSetOUT(parent1ID,8) + segmentSetOUT(parent2ID,8) ~= ...
                            segmentSetOUT(childID,8) ) && (~any(segmentErr{1} == parent1ID) && ...
                            ~any(segmentErr{1} == parent2ID)) )
                        %There is an intensity size discrepancy.
                        segmentErr{1}(sum(~isnan(segmentErr{1}))+1) = parent1ID;
                    end              
                end %if isnan child 1 of parent 1

                segmentSetOUT(childID,9) = segmentSetOUT(parent1ID,9) + ...
                    segmentSetOUT(parent2ID,9);  
                callCount = callCount + 1;
                [segmentSetOUT,segmentErr,callCount] = processSegmentEvents(nextEventTime,childID,segmentSetOUT,...
                    segmentErr,seqOfEvents,tracksAmp,intensityInfo,callCount);
                callCount = callCount - 1;

                %Check amplitude from events after call - not needed but
                %keep anyway
                %{
                if (segmentSetOUT(parent1ID,9) + segmentSetOUT(parent2ID,9) ~= ...
                        segmentSetOUT(childID,9) && (~any(segmentErr{2} == parent1ID) && ...
                            ~any(segmentErr{2} == parent2ID)) )
                    %There is an event size discrepancy that can only be
                    %corrected using intensity values.  This is  due to a 
                    %merge being followed by a split - in previous
                    %code, this is the case handled in the split trace-back 
                    %block where parent originated from a merge event
                    segmentErr{2}(sum(~isnan(segmentErr{2}))+1) = parent1ID;
                end                                        
                  %}  
            end %If nextEventType
            
        end %if nextEventRow not empty
        
    %end %if segmentSetIN has NaN fields
    
end %function

