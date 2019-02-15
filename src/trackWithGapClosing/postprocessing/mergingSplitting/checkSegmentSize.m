function ampMatch = checkSegmentSize(segID,segmentSetIN,direction,sizeType,checkAll)
%CHECKSEGMENTSIZE performs segment size checks for events involving parent
%or children segments of the given segment.
%
%   INPUT:
%           segID:          id for the reference segment
%           segmetSetIN:    table of segments and events constructed by
%                           processSegmentEvents
%           direction:      if a value of 1 given, then events for children
%                           of segID will be checked and if value of -1
%                           given, events for parents of segID will be
%                           checked.
%           sizeType:       this function is used by
%                           repairSizeFromIntensity and 
%                           repairSizeFromEvents for the two types of sizes
%                           collected. Possible values are 1 or 2, for
%                           intensity or event sizes, respectively.
%           checkAll:       a boolean value indicating whther the size
%                           check should propagate in the specified
%                           direction until the start or end is reached.
%   OUTPUT:
%           ampMatch:       result of the size check:
%                           NaN, if forward/backward events do not exist 
%                           for segID.
%                           0, if forward/backward events exist but there 
%                           is a size problem
%                           1, if forward/backward events exist and there 
%                           is no size problem
%
%   Robel Yirdaw, December 2013
%

    %Initialize return value - default is NaN
    ampMatch = NaN;
    %Incoming segment is parent
    parentID = segID;
    %Determine which amplitude is checked. 1 = intensity, 2 = event.
    sizeIndx = NaN;
    if (sizeType == 1)
        sizeIndx = 8;
    elseif (sizeType == 2)
        sizeIndx = 9;
    end
    %The following two are for checkAll = 1
    outStr = '';
    nextSegID = [NaN;NaN];
    
    parentSize = NaN;
    if (direction == -1)
        %Backwards check
        if (~isnan(segmentSetIN(parentID,1)) && isnan(segmentSetIN(parentID,2)) && ...
                ~isnan(segmentSetIN(parentID,5)))
            %Parent was created via split - either it or its sibling must
            %have size 1
            if ( ( (segmentSetIN(parentID,sizeIndx) == 1) || ...
                    (segmentSetIN(segmentSetIN(parentID,5),sizeIndx) == 1) ) && ...
                    ( (segmentSetIN(segmentSetIN(parentID,1),sizeIndx) > 0) && ...
                      (segmentSetIN(segmentSetIN(parentID,5),sizeIndx) > 0) )  )
                parentSize = segmentSetIN(segmentSetIN(parentID,1),sizeIndx) - ...
                    segmentSetIN(segmentSetIN(parentID,5),sizeIndx);
            else
                parentSize = 0;
            end
            if(checkAll)
                outStr = sprintf('\n(%d:%d) -> (%d:%d) + (%d:%d). ',segmentSetIN(parentID,1),...
                    segmentSetIN(segmentSetIN(parentID,1),sizeIndx),parentID,segmentSetIN(parentID,sizeIndx),...
                    segmentSetIN(parentID,5),segmentSetIN(segmentSetIN(parentID,5),sizeIndx));
                %will continue with parent's parent
                nextSegID(1) = segmentSetIN(parentID,1);
            end
        elseif (~isnan(segmentSetIN(parentID,1)) && ~isnan(segmentSetIN(parentID,2)))
            %Parent was created via merge - one of the parent's parents
            %must have size 1
            if ( ( (segmentSetIN(segmentSetIN(parentID,1),sizeIndx) == 1) || ...
                    (segmentSetIN(segmentSetIN(parentID,2),sizeIndx) == 1) ) &&...
                    ( (segmentSetIN(segmentSetIN(parentID,1),sizeIndx) > 0) && ...
                      (segmentSetIN(segmentSetIN(parentID,2),sizeIndx) > 0) ) )
                parentSize = segmentSetIN(segmentSetIN(parentID,1),sizeIndx) + ...
                    segmentSetIN(segmentSetIN(parentID,2),sizeIndx);
            else
                parentSize = 0;
            end
            if(checkAll)
                outStr = sprintf('\n(%d:%d) + (%d:%d) -> (%d:%d). ',segmentSetIN(parentID,1),...
                    segmentSetIN(segmentSetIN(parentID,1),sizeIndx),segmentSetIN(parentID,2),...
                    segmentSetIN(segmentSetIN(parentID,2),sizeIndx),parentID,segmentSetIN(parentID,sizeIndx));
                %will continue with both parents of parent
                nextSegID(1) = segmentSetIN(parentID,1);
                %{
                if (segmentSetIN(segmentSetIN(parentID,1),1) ~= ...
                        segmentSetIN(segmentSetIN(parentID,2),1) )
                    nextSegID(2) = segmentSetIN(parentID,2);
                end
                %}
            end            
        end
    elseif (direction == 1)
        %Forward check
        if (~isnan(segmentSetIN(parentID,3)) && ~isnan(segmentSetIN(parentID,4)))
            %Parent is splitting - one of the children must be 1
            if ( ((segmentSetIN(segmentSetIN(parentID,3),sizeIndx) == 1) || ...
                    (segmentSetIN(segmentSetIN(parentID,4),sizeIndx) == 1) ) && ...
                    ( (segmentSetIN(segmentSetIN(parentID,3),sizeIndx) > 0) && ...
                      (segmentSetIN(segmentSetIN(parentID,4),sizeIndx) > 0) ) )
                parentSize = segmentSetIN(segmentSetIN(parentID,3),sizeIndx) + ...
                    segmentSetIN(segmentSetIN(parentID,4),sizeIndx);
            else
                parentSize = 0;
            end
            if (checkAll)
                outStr = sprintf('\n(%d:%d) -> (%d:%d) + (%d:%d). ',parentID,segmentSetIN(parentID,sizeIndx),...
                    segmentSetIN(parentID,3),segmentSetIN(segmentSetIN(parentID,3),sizeIndx),...
                    segmentSetIN(parentID,4),segmentSetIN(segmentSetIN(parentID,4),sizeIndx));
                %will continue with both children
                nextSegID(1) = segmentSetIN(parentID,3);
                %{
                if (segmentSetIN(segmentSetIN(parentID,3),3) ~= ...
                        segmentSetIN(segmentSetIN(parentID,4),3) )
                    nextSegID(2) = segmentSetIN(parentID,4);
                end
                %}
            end
        elseif (~isnan(segmentSetIN(parentID,3)) && isnan(segmentSetIN(parentID,4)) && ...
                ~isnan(segmentSetIN(parentID,6)))
            %Parent is merging - either it or its partner must have size 1
            if ( ( (segmentSetIN(parentID,sizeIndx) == 1) || ...
                    (segmentSetIN(segmentSetIN(parentID,6),sizeIndx) == 1) ) && ...
                    ( (segmentSetIN(segmentSetIN(parentID,3),sizeIndx) > 0) && ...
                      (segmentSetIN(segmentSetIN(parentID,6),sizeIndx) > 0) ) )
                parentSize = segmentSetIN(segmentSetIN(parentID,3),sizeIndx) - ...
                    segmentSetIN(segmentSetIN(parentID,6),sizeIndx);
            else
                parentSize = 0;
            end
            if (checkAll)
                outStr = sprintf('\n(%d:%d) + (%d:%d) -> (%d:%d). ',parentID,segmentSetIN(parentID,sizeIndx),...
                    segmentSetIN(parentID,6),segmentSetIN(segmentSetIN(parentID,6),sizeIndx),...
                    segmentSetIN(parentID,3),segmentSetIN(segmentSetIN(parentID,3),sizeIndx));
                %will continue with child
                nextSegID(1) = segmentSetIN(parentID,3);
            end            
        end
    end
    
    %if parentSize has not been computed, it will be NaN. This can happen
    %when the parent segment doesn't have parents or children.
    if (~isnan(parentSize))
        if ( (segmentSetIN(parentID,sizeIndx) > 0) && (parentSize > 0) )
            ampMatch = (segmentSetIN(parentID,sizeIndx) == parentSize);
        else
            ampMatch = 0;
        end
    end
    
    %fprintf(outStr);
    if (checkAll)
        if (ampMatch ~= 1)
            fprintf(outStr);
        end
        if (~isnan(nextSegID(1)))
            checkSegmentSize(nextSegID(1),segmentSetIN,direction,sizeType,checkAll);
            %fprintf('\n Return on %d. ',nextSegID(1));
            if (~isnan(nextSegID(2)))
                checkSegmentSize(nextSegID(2),segmentSetIN,direction,sizeType,checkAll);
                %fprintf('  Return on %d. ',nextSegID(2));
            end
        end
    end
    
end %function
