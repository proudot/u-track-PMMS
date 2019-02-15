function [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(segID,segmentSetIN,direction,notRepEveSizeErr)
%REPAIRSIZEFROMEVENTS corrects segment sizes that were determined 
%exclusively from the sequence of merging and splitting events. Only on
%end points (beginning or ending of events) will intensity amplitude be
%used to assign sizes for segments. This function performs repairs
%recursively.
%
%   INPUT:
%           segID:              ID of segment currently being processed
%           segmentSetIN:       table with each row corresponding to a
%                               segment and the following columns as
%                               constructed by processSegmentEvents:
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
%           direction:          if a value of 1 given, then the repair will
%                               progress forwards with the childrend of segID
%                               segID and if a value of -1 is given, then 
%                               the repair continues backwards, with the 
%                               parents of segID.
%           notRepEveSizeErr:   a 3x1 cell containing the following
%                               elements for those segments whose repair
%                               failed:
%                               1) The segment ID
%                               2) direction the repair was moving
%                               3) a string stating the error
%
%   OUTPUT:
%           segementSetOUT:     updated table as detailed above
%           notRepEveSizeErr:   final 3x1 cell as detailed above
%
%   NOTE: this function was never used and as a result is not fully
%   developed and tested. 
%
%   Robel Yirdaw, November 2013.
%
    
    segmentSetOUT = segmentSetIN;

    if (direction == 0)
        %This block here sets up the sequence of corrections that will be
        %needed for the incoming segment. It also works to prevent the case
        %where a segment is in the error list but has been corrected via
        %another segment's correction earlier. In this case, the coditions
        %checked below will prevent the directional search to be initiated.
        
        %Starting segment
        parentID = segID;
        if (~isnan(segmentSetIN(parentID,3)) && ~isnan(segmentSetIN(parentID,4)))    
            %This is a split event. Okay to continue.
            child1Size = segmentSetIN(segmentSetIN(parentID,3),9);
            child2Size = segmentSetIN(segmentSetIN(parentID,4),9);
            parentSize = segmentSetIN(parentID,9);
            
            if (parentSize < child1Size + child2Size)                
                %Update parent's amplitude
                segmentSetOUT(parentID,9) = child1Size + child2Size;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %START BACKWARDS CORRECTIONS                        
                [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(parentID,segmentSetOUT,-1,notRepEveSizeErr);                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif (parentSize > child1Size + child2Size)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %START FORWARDS CORRECTIONS
                [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(parentID,segmentSetIN,1,notRepEveSizeErr);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end %if parent's amplitude < childrens'

        end %incoming segment is splitting
        
    elseif (direction == -1)
        %Traversing backwards - incoming segment is child segment.
        childID = segID;
        parent1ID = segmentSetIN(childID,1);
        parent2ID = segmentSetIN(childID,2);
        %Could have been created via merge or split
        
        if ((~isnan(parent1ID) && isnan(parent2ID)) && (segmentSetIN(parent1ID,9) ~= ...
                segmentSetIN(childID,9) + segmentSetIN(segmentSetIN(childID,5),9) ) )
            %The incoming child segment created via split; has sibling
            siblingID = segmentSetIN(childID,5);
            
            %Need to prevent splitting of clusters
            if ((segmentSetIN(childID,9) > 1) && segmentSetIN(siblingID,9) > 1)
                %The update on child amplitude does not work. Since
                %Revert. 
                segmentSetOUT(childID,9) = segmentSetOUT(parent1ID,9) - ...
                    segmentSetOUT(siblingID,9);
            else            
                %Update the parent's amplitude            
                segmentSetOUT(parent1ID,9) = segmentSetIN(childID,9) + ...
                    segmentSetIN(siblingID,9);            
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Continue upwards with parent, if it has one or two parents.
                %Otherwise, this is an initial segment that starts with a
                %split.
                if (~isnan(segmentSetOUT(parent1ID,1)))
                    [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(parent1ID,segmentSetOUT,-1,notRepEveSizeErr);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
        elseif ((~isnan(parent1ID) && ~isnan(parent2ID)) && (segmentSetIN(childID,9) ~= ...
                    segmentSetIN(parent1ID,9) + segmentSetIN(parent2ID,9)) )
            %The incoming child segment created via merge
            
            if (segmentSetIN(childID,9) == 1)
                %The update on child amplitude does not work. Revert. A
                %merge can not produce a segment of size 1.
                segmentSetOUT(childID,9) = segmentSetIN(parent1ID,9) + ...
                    segmentSetIN(parent2ID,9);
            elseif (segmentSetIN(childID,9) == 2)
                segmentSetOUT(parent1ID,9) = 1;
                segmentSetOUT(parent2ID,9) = 1;
            else
                if ((segmentSetIN(parent1ID,9) == 1) && (segmentSetIN(parent2ID,9) == 1))
                    %The parents have the same size.  Use intensity.  
                    %{
                    fprintf('\n(%d,%d) Merging segments have equal size: (%d:%d) + (%d:%d) -> (%d:%d). ',...
                        segID,direction,parent1ID,segmentSetIN(parent1ID,9),parent2ID,segmentSetIN(parent2ID,9),...
                        childID,segmentSetIN(childID,9));
                    %}
                    parent12ID = [parent1ID;parent2ID];                    
                    parent12Amp = [segmentSetIN(parent1ID,8);segmentSetIN(parent2ID,8)];
                    %identify which segment has higher intensity amplitude,
                    segmentBrighter = find(parent12Amp == max(parent12Amp));
                    if (length(segmentBrighter) ~= 1)
                        outStr = sprintf('\n(%d,%d) Unable to correct event sizes - merging segments (%d & %d) have the same intensity amplitudes. ',...
                            segID,direction,parent1ID,parent2ID);
                        tempIndx = sum(~isnan(notRepEveSizeErr{1})) + 1;
                        notRepEveSizeErr{1}(tempIndx) = childID;
                        notRepEveSizeErr{2}(tempIndx) = direction;
                        notRepEveSizeErr{3}{tempIndx} = outStr;
                        fprintf(outStr);
                    else
                        %{
                        fprintf('\n(%d,%d) Adjusting *merge* event amplitude using intensity amplitude:  (%d:%d) + (%d:%d) -> (%d:%d).  ',...
                            segID,direction,parent1ID,segmentSetIN(parent1ID,8),parent2ID,segmentSetIN(parent2ID,8),...
                            childID,segmentSetIN(childID,8));
                        %}
                        %segmentBrighter = segmentBrighter(1);
                        segmentDimmer = setdiff([1;2],segmentBrighter);
                        parentBrighterID = parent12ID(segmentBrighter);
                        parentDimmerID = parent12ID(segmentDimmer);

                        %update the amplitude values for the merging parents of
                        %the current parent
                        %segmentSetOUT(parentDimmerID,9) = parent12Amp(segmentDimmer);
                        %segmentSetOUT(parentBrighterID,9) = parent12Amp(segmentBrighter);
                        %Adjust by subtracting 1 because clusters do not merge.
                        segmentSetOUT(parentBrighterID,9) = segmentSetIN(childID,9) - 1;
                        segmentSetOUT(parentDimmerID,9) = 1;
                    end
                elseif (segmentSetIN(parent1ID,9) > 1)
                    segmentSetOUT(parent1ID,9) = segmentSetIN(childID,9) - 1;
                    segmentSetOUT(parent2ID,9) = 1;
                elseif (segmentSetIN(parent2ID,9) > 1)
                    segmentSetOUT(parent2ID,9) = segmentSetIN(childID,9) - 1;
                    segmentSetOUT(parent1ID,9) = 1;               
                end
            end %if incoming child segment size

            %Finished settings sizes. Continue as possible.
            if (~isnan(segmentSetIN(parent1ID,1)) && (segmentSetOUT(parent1ID,9) ~= ...
                    segmentSetIN(parent1ID,9)) )
                %Parent 1 continues backwards
                [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(parent1ID,segmentSetOUT,-1,notRepEveSizeErr);
                %Check the returned values
                if ( (segmentSetOUT(childID,9) ~= segmentSetOUT(parent1ID,9) + ...
                        segmentSetOUT(parent2ID,9)) && (segmentSetIN(parent1ID,9) ~= ...
                        segmentSetOUT(parent2ID,9)) )
                    %Update did not work. Try reversing values. 
                    %111113
                    %NOTE: unlike intensity sizes, here there is no need to
                    %check whether parent2 has parent's and if the swap
                    %affects parent2's parents. If parent2 has parents, it
                    %will be traversed below. If not, then the swap is ok.
                    tempParent1Size = segmentSetOUT(parent1ID,9);
                    segmentSetOUT(parent1ID,9) = segmentSetOUT(parent2ID,9);
                    segmentSetOUT(parent2ID,9) = tempParent1Size;
                    %Call with reversed values and parent1 ID
                    [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(parent1ID,segmentSetOUT,-1,notRepEveSizeErr);
                end
            end
            if (~isnan(segmentSetIN(parent2ID,1)) && ...
                    (segmentSetOUT(parent2ID,9) ~= segmentSetIN(parent2ID,9)) && ...
                    ( (segmentSetIN(parent1ID,1) ~= segmentSetIN(parent2ID,1)) || ...
                      ( (segmentSetIN(parent1ID,1) == segmentSetIN(parent2ID,1)) && ...
                        (segmentSetOUT(parent1ID,9) == segmentSetIN(parent1ID,9)) ) ) )
                %Parent 2 continues backwards and has a different
                %parent than parent 1
                [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(parent2ID,segmentSetOUT,-1,notRepEveSizeErr);
                %Check the returned values
                if ((segmentSetOUT(childID,9) ~= segmentSetOUT(parent1ID,9) + ...
                        segmentSetOUT(parent2ID,9)) && isnan(segmentSetIN(parent1ID,1)) && ...
                        (segmentSetOUT(parent1ID,9) ~= segmentSetIN(parent2ID,9)) )
                    %Update did not work. Try reversing values.
                    tempParent2Size = segmentSetOUT(parent2ID,9);                        
                    segmentSetOUT(parent2ID,9) = segmentSetOUT(parent1ID,9);
                    segmentSetOUT(parent1ID,9) = tempParent2Size;
                    %Call with reversed values and parent2 ID
                    segmentSetOUT = repairSizeFromEvents(parent2ID,segmentSetOUT,-1,notRepEveSizeErr);
                end
            end
        end %if incoming seg from split or merge
        
    elseif (direction == 1)
        %Traversing forwards - incoming segment is parent segment.
        parentID = segID;
        %Parent can be splitting or merging
        child1ID = segmentSetIN(parentID,3);
        child2ID = segmentSetIN(parentID,4);
        if ((~isnan(child1ID) && isnan(child2ID) && ~isnan(segmentSetIN(parentID,6)) ) && ...
                (segmentSetIN(child1ID,9) ~= segmentSetIN(parentID,9) + ...
                    segmentSetIN(segmentSetIN(parentID,6),9) ) )
            %Parent is merging - has only one child and has partner.
            partnerID = segmentSetIN(parentID,6);
            %Check that clusters are not merging - one has to have size 1.
            if ((segmentSetIN(parentID,9) > 1) && (segmentSetIN(partnerID,9) > 1))
                %The update on the incoming segment did not work. Revert.
                segmentSetOUT(parentID,9) = segmentSetIN(child1ID,9) - ...
                    segmentSetIN(partnerID,9);
            else
                segmentSetOUT(child1ID,9) = segmentSetIN(parentID,9) + ...
                    segmentSetIN(partnerID,9);
                %Continue with child if child has at least one child
                if (~isnan(segmentSetIN(child1ID,3)))
                    [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(child1ID,segmentSetOUT,1,notRepEveSizeErr);
                end
            end %If incoming updated segmet ok
        
        elseif ( (~isnan(child1ID) && ~isnan(child2ID)) && (segmentSetIN(child1ID,9) + ...
                    segmentSetIN(child2ID,9) ~= segmentSetIN(parentID,9)) )
            %Parent is splitting - has two children
            %If the incoming parent segment has size 1, can not continue,
            %since it is splitting.
            if (segmentSetIN(parentID,9) == 1)
                %The update on the incoming segment did not work. Revert.
                segmentSetOUT(parentID,9) = segmentSetIN(child1ID,9) + ...
                    segmentSetIN(child2ID,9);
            elseif (segmentSetIN(parentID,9) == 2)
                segmentSetOUT(child1ID,9) = 1;
                segmentSetOUT(child2ID,9) = 1;
            else
                %parent's size > 2.
                child1IntAmp = segmentSetIN(child1ID,8);
                child2IntAmp = segmentSetIN(child2ID,8);   
                parentSize = segmentSetIN(parentID,9);
                if ((segmentSetIN(child1ID,9) == 1) && (segmentSetIN(child2ID,9) == 1))
                    %The splitting segments have the same size. This
                    %can only happen if they both have size = 1. Thus,
                    %impossible to pick one path to modify simply based
                    %on sequence of events.   
                    %{
                    fprintf('\n(%d,%d) Splitting segments have equal size: (%d:%d) -> (%d:%d) + (%d:%d).  ',...
                        segID,direction,parentID,parentSize,child1ID,segmentSetIN(child1ID,9),...
                        child2ID,segmentSetIN(child2ID,9));
                    %}
                    %Use intensity amplitude to pick the brighter
                    %segment to continue with
                    child12ID = [child1ID;child2ID];
                    child12IntAmp = [child1IntAmp;child2IntAmp];
                    %identify which segment has lower amplitude and which
                    %has higher amplitude                        
                    segmentDimmer = find(child12IntAmp == min(child12IntAmp));
                    if (length(segmentDimmer) ~= 1)
                        outStr = sprintf('\n(%d,%d) Unable to correct event sizes - segments (%d & %d) have the same intensity amplitudes. ',...
                            segID,direction,child1ID,child2ID);
                        tempIndx = sum(~isnan(notRepEveSizeErr{1})) + 1;
                        notRepEveSizeErr{1}(tempIndx) = parentID;
                        notRepEveSizeErr{2}(tempIndx) = direction;
                        notRepEveSizeErr{3}{tempIndx} = outStr;
                        fprintf(outStr);                        
                    else
                        %{
                        fprintf('\n(%d,%d) Adjusting *split* event amplitude using intensity amplitude: (%d:%d) -> (%d:%d) + (%d:%d). ',...
                            segID,direction,parentID,segmentSetIN(segID,8),child1ID,segmentSetIN(child1ID,8),...
                        child2ID,segmentSetIN(child2ID,8));      
                        %}
                        segmentDimmer = segmentDimmer(1);
                        segmentBrighter = setdiff([1;2],segmentDimmer);      

                        segmentSetOUT(child12ID(segmentBrighter),9) = parentSize - 1;
                        segmentSetOUT(child12ID(segmentDimmer),9) = 1;
                    end  
                elseif ((segmentSetIN(child1ID,9) > 1))
                    segmentSetOUT(child1ID,9) = parentSize - 1;   
                    segmentSetOUT(child2ID,9) = 1;
                elseif ((segmentSetIN(child2ID,9) > 1))
                    segmentSetOUT(child2ID,9) = parentSize - 1;                                                        
                    segmentSetOUT(child1ID,9) = 1;                
                end %if childrens sizes
            end %if parent's size

            %Done setting up sizes
            if (~isnan(segmentSetIN(child1ID,3)) && (segmentSetOUT(child1ID,9) ~= ...
                    (segmentSetIN(child1ID,9)) )  )
                 [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(child1ID,segmentSetOUT,1,notRepEveSizeErr);
                %Check the returned values
                if ( (segmentSetOUT(parentID,9) ~= segmentSetOUT(child1ID,9) + ...
                        segmentSetOUT(child2ID,9)) && (segmentSetIN(child1ID,9) ~= ...
                    (segmentSetOUT(child2ID,9)) )  ) 
                    %Update did not work. Try reversing values.
                    tempChild1Size = segmentSetOUT(child1ID,9);                        
                    segmentSetOUT(child1ID,9) = segmentSetOUT(child2ID,9);
                    segmentSetOUT(child2ID,9) = tempChild1Size;
                    %Call with reversed values and child1 ID
                    [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(child1ID,segmentSetOUT,1,notRepEveSizeErr);
                end

            end
            if (~isnan(segmentSetIN(child2ID,3)) && ...
                    (segmentSetOUT(child2ID,9) ~= segmentSetIN(child2ID,9)) && ...
                    ( (segmentSetIN(child1ID,3) ~= segmentSetIN(child2ID,3)) || ...
                        ( (segmentSetIN(child1ID,3) == segmentSetIN(child2ID,3)) && ...
                          (segmentSetOUT(child1ID,9) == segmentSetIN(child1ID,9)) )  ) )
                [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(child2ID,segmentSetOUT,1,notRepEveSizeErr);
                if ((segmentSetOUT(parentID,9) ~= segmentSetOUT(child1ID,9) + ...
                        segmentSetOUT(child2ID,9)) && (isnan(segmentSetIN(child1ID,3))) && ...
                        (segmentSetOUT(child1ID,9) ~= segmentSetIN(child2ID,9)) )
                    %Update did not work. Try reversing values.
                    tempChild2Size = segmentSetOUT(child2ID,9);                        
                    segmentSetOUT(child2ID,9) = segmentSetOUT(child1ID,9);
                    segmentSetOUT(child1ID,9) = tempChild2Size;
                    %Call with reversed values and child1 ID
                    [segmentSetOUT,notRepEveSizeErr] = repairSizeFromEvents(child2ID,segmentSetOUT,1,notRepEveSizeErr);
                end
            end
        
        end %if parent (incoming segment) merging or splitting
        
    end %if direction

end %function

