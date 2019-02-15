function [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(...
    segID,segmentSetIN,checkAll,notRepIntSizeErr,intensityInfo)
%REPAIRSIZEFROMINTENSITY corrects errors in segment sizes that were 
%determined from intensity amplitudes. This function performs repairs 
%recursively.
%
%   Repair of the incorrect size of a segment begins by first looking at
%   the event the segment participates in, if any. Depending on the type of
%   event, split or merge, the child segment(s) will be checked for size.
%   These themselves could have size errors. In that case, the sizes for
%   each child segment will be checked, and so on, until a segment with
%   a correct size assignment is found. This is the recursive part of the
%   function. Each call will then return repairing the child sizes until
%   the very first segment is reached and repaired.
%
%   In the event that sizes can not be corrected by looking at the
%   respetive events and sizes for the other segmetns involved, sizes will
%   be determined probablistically. Thus, there should be no case where a
%   size error on child segments returns unrepaired. If the probabilistic
%   approach is not used, recursive repairs on child segments can return
%   without repairing the segments. In that case, the separate cases listed
%   in the function need to be able to handle these scenarios. The
%   development version of this code includes all of those scenarios for
%   each split and merge case type.
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
%           checkAll:           a boolean value indicating whether repairs
%                               should continue until the end of events
%           notRepIntSizeErr:   a 5x1 cell initialized in 
%                               aggregStateFromCompTracks_new containing 
%                               the following elements for those segments
%                               whose sizes could not be repaired:
%                               1) ID of segment not repaired
%                               2) a cell with a string indicating error
%                               case type
%                               3) a cell with a string stating the repair
%                               outcome
%                               4) a cell with a string for the event type,
%                               either a split or merge
%                               5) the amount of discrepancy, calculated
%                               based on the type of event
%           intensityInfo:      intensity quantum as defined in main sim.
%
%   OUTPUT:
%           segementSetOUT:     updated table as detailed above
%           notRepIntSizeErr:   final 3x1 cell as detailed above
%
%
%   Robel Yirdaw, November 2013
%       Modified, December 2014
%   

    %Set to display messages to screen
    verbose = true;
    %Copy incoming segmentSet
    segmentSetOUT = segmentSetIN;
    
    %Begin
    if ( isnan(segmentSetIN(segID,3)) && isnan(segmentSetIN(segID,4)) )
        fprintf('\nSegment %d has no children. ',segID);
    elseif ( ~isnan(segmentSetIN(segID,3)) && isnan(segmentSetIN(segID,4)) && ...
            isnan(segmentSetIN(segID,6)) )
        fprintf('\nError: segment %d merging with unknown segment. ',segID);
    elseif ( ~isnan(segmentSetIN(segID,3)) && ~isnan(segmentSetIN(segID,4)) )
        %Split event.
        parentID = segID;
        child1ID = segmentSetIN(parentID,3);
        child2ID = segmentSetIN(parentID,4);
        parentSize = segmentSetIN(parentID,8);
        child1Size = segmentSetIN(child1ID,8);
        child2Size = segmentSetIN(child2ID,8);

        if ( (child1Size + child2Size ~= parentSize) ||... 
             (child1Size == 0) || (child2Size == 0) || (parentSize == 0) ||...
                 ((child1Size > 1) && (child2Size > 1))  )
             %{
            fprintf('\n(%d:%d) -> (%d:%d) + (%d:%d). ',parentID,parentSize,...
                child1ID,child1Size,child2ID,child2Size);
             %}
            %Try correcting child1
            [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(child1ID,segmentSetIN,0,notRepIntSizeErr,intensityInfo);
            child1Size_rev = segmentSetOUT(child1ID,8);
            
            if ( (child1Size_rev + child2Size ~= parentSize) ||... 
                 (child1Size_rev == 0) || (child2Size == 0) || (parentSize == 0) ||...
                 ((child1Size_rev > 1) && (child2Size > 1)) )
                %Try correcting child2. Traverse only if the two child
                %segments do not merge back with eachother. (Rearrange).
                if ( segmentSetIN(child2ID,3)~= segmentSetIN(child1ID,3))
                    [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(child2ID,segmentSetOUT,0,notRepIntSizeErr,intensityInfo);
                end
                child2Size_rev = segmentSetOUT(child2ID,8);
                
                if ( (child1Size + child2Size_rev == parentSize) && (child1Size > 0) && ...
                (child2Size_rev > 0) && (parentSize > 0) &&...
                ((child1Size == 1) || (child2Size_rev == 1)) )                    
                    fprintf('Repaired (s1): (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,parentSize,...
                        child1ID,child1Size,child2ID,child2Size_rev);                                        
                elseif ( (child1Size_rev + child2Size_rev == parentSize) && (child1Size_rev > 0) && ...
                (child2Size_rev > 0) && (parentSize > 0) &&...
                ((child1Size_rev == 1) || (child2Size_rev == 1))  )
                    segmentSetOUT(segmentSetOUT(segID,3),8) = child1Size_rev;                    
                    fprintf('Repaired (s2): (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,parentSize,...
                        child1ID,child1Size_rev,child2ID,child2Size_rev);                    
                else
                    %Both child segments checked. Try to repair via swap.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    swapFxd = 0;
                    if ( (child1Size_rev > 0) && (child2Size_rev > 0) && (parentSize > 0) )
                        
                        parentSibSize = NaN;
                        child1PartSize = NaN;
                        child2PartSize = NaN;                    
                        if (~isnan(segmentSetOUT(parentID,5)) && isnan(segmentSetOUT(segmentSetOUT(parentID,5),3)) &&...
                                (segmentSetOUT(segmentSetOUT(parentID,5),8) > 0) )
                            parentSibSize = segmentSetOUT(segmentSetOUT(parentID,5),8);
                        end
                        if (~isnan(segmentSetOUT(child1ID,6)) && isnan(segmentSetOUT(segmentSetOUT(child1ID,6),1)) &&...
                                (segmentSetOUT(segmentSetOUT(child1ID,6),8) > 0) )
                            child1PartSize = segmentSetOUT(segmentSetOUT(child1ID,6),8);
                        end
                        if (~isnan(segmentSetOUT(child2ID,6)) && isnan(segmentSetOUT(segmentSetOUT(child2ID,6),1)) && ...
                                (segmentSetOUT(segmentSetOUT(child2ID,6),8) > 0) )
                            child2PartSize = segmentSetOUT(segmentSetOUT(child2ID,6),8);
                        end
                        child1Vec = [child1Size_rev; child1PartSize; child1Size_rev; child1PartSize];
                        child2Vec = [child2Size_rev; child2Size_rev; child2PartSize; child2PartSize];
                        childSumVec = child1Vec + child2Vec;
                        correctSum = find(childSumVec == parentSibSize,1,'first');
                        if (~isempty(correctSum))
                            %A swap that works found. Assign values.
                            %Swap parent's size with it's sibling's
                            segmentSetOUT(parentID,8) = parentSibSize;
                            segmentSetOUT(segmentSetOUT(parentID,5),8) = parentSize;
                            %Swap child1's size with it's partner's
                            if (child1Vec(correctSum) == child1PartSize)
                                segmentSetOUT(segmentSetOUT(child1ID,6),8) = child1Size_rev;
                                segmentSetOUT(child1ID,8) = child1Vec(correctSum);
                            end
                            %Swap child2's size with it's partner's
                            if (child2Vec(correctSum) == child2PartSize)
                                segmentSetOUT(segmentSetOUT(child2ID,6),8) = child2Size_rev;
                                segmentSetOUT(child2ID,8) = child2Vec(correctSum);
                            end
                            swapFxd = 1;    
                            fprintf('Repaired (swap sa): (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,segmentSetOUT(parentID,8),...
                                child1ID,segmentSetOUT(child1ID,8),child2ID,segmentSetOUT(child2ID,8));                    
                        else
                            correctSum = find(childSumVec(2:4) == parentSize,1,'first');
                            if (~isempty(correctSum))
                                %A swap that works found. Assign values.
                                %Only swapping children sizes here.
                                %Swap child1's size with it's partner's
                                if (child1Vec(correctSum) == child1PartSize)
                                    segmentSetOUT(segmentSetOUT(child1ID,6),8) = child1Size_rev;
                                    segmentSetOUT(child1ID,8) = child1Vec(correctSum);
                                end
                                %Swap child2's size with it's partner's
                                if (child2Vec(correctSum) == child2PartSize)
                                    segmentSetOUT(segmentSetOUT(child2ID,6),8) = child2Size_rev;
                                    segmentSetOUT(child2ID,8) = child2Vec(correctSum);
                                end
                                swapFxd = 1;    
                                fprintf('Repaired (swap sb): (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,segmentSetOUT(parentID,8),...
                                    child1ID,segmentSetOUT(child1ID,8),child2ID,segmentSetOUT(child2ID,8));                    
                            end
                        end
                    end %if child and parent sizes > 0
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if (~swapFxd)
                        %Attempt to repair via assignment
                        %First, check if parent's prior events are ok. 
                        parentChk = checkSegmentSize(parentID,segmentSetOUT,-1,1,0);
                        if (isnan(parentChk) && isnan(segmentSetOUT(child1ID,3)) && ...
                                isnan(segmentSetOUT(child2ID,3)) )
                            %case 1
                            
                            segmentSetOUT = assignSizeFromProb([parentID;NaN],parentSize,...
                                [child1ID;child2ID],segmentSetIN,segmentSetOUT,intensityInfo);      
                            
                            fprintf('Repaired (1): (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,segmentSetOUT(parentID,8),...
                                child1ID,segmentSetOUT(child1ID,8),child2ID,segmentSetOUT(child2ID,8));                                     

                        elseif (isnan(parentChk) && isnan(segmentSetOUT(child1ID,3)) && ...
                                ~isnan(segmentSetOUT(child2ID,3)) )
                            %case 4
                            tempChild2Chk = checkSegmentSize(child2ID,segmentSetOUT,1,1,0);
                            tempRepFlag = 0;
                            if (tempChild2Chk == 1)
                                %child ok.
                                if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                    segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                    tempRepFlag = 1;
                                elseif ((child1Size_rev > 1) && (child2Size_rev > 1)) 
                                    segmentSetOUT(child1ID,8) = 1;
                                    segmentSetOUT(parentID,8) = child2Size_rev + 1;
                                    tempRepFlag = 1;
                                end
                            end %if tempChild2Chk failed
                            
                            clear tempChild2Chk
                            
                            notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'4',verbose);                            
                            
                        elseif (isnan(parentChk) && ~isnan(segmentSetOUT(child1ID,3)) && ...
                                isnan(segmentSetOUT(child2ID,3)) )
                            %case 7
                            tempChild1Chk = checkSegmentSize(child1ID,segmentSetOUT,1,1,0);
                            tempRepFlag = 0;
                            if (tempChild1Chk == 1)
                                %child ok.
                                if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                    segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                    tempRepFlag = 1;
                                elseif ((child1Size_rev > 1) && (child2Size_rev > 1)) 
                                    segmentSetOUT(child2ID,8) = 1;
                                    segmentSetOUT(parentID,8) = child1Size_rev + 1;
                                    tempRepFlag = 1;
                                end
                                
                            end %if tempChild1Chk failed
                            
                            clear tempChild1Chk
                            
                            notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'7',verbose);                             
                            
                        elseif (isnan(parentChk) && ~isnan(segmentSetOUT(child1ID,3)) && ...
                                ~isnan(segmentSetOUT(child2ID,3)) )
                            %case 10
                            tempChild1Chk = checkSegmentSize(child1ID,segmentSetOUT,1,1,0);
                            tempChild2Chk = checkSegmentSize(child2ID,segmentSetOUT,1,1,0);
                            tempRepFlag = 0;    
                            
                            if ( (tempChild1Chk == 1) && (tempChild2Chk == 1) )
                                if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                    segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                    tempRepFlag = 1;                   
                                else
                                    %Both > 1. The lower one must be set to
                                    %1 and then repaired. Can't change the
                                    %other one because it's check = 1.
                                    childIDs = [child1ID;child2ID];
                                    %childSizes = [child1Size_rev;child2Size_rev];
                                    childInts = [segmentSetIN(child1ID,7);segmentSetIN(child2ID,7)];
                                    %tempBigger = find((childInts == max(childInts)));
                                    smallerChildID = childIDs((childInts == min(childInts)));
                                    biggerChildID = childIDs(childIDs ~= smallerChildID);
                                    segmentSetOUT(smallerChildID(1),8) = 1;
                                    segmentSetOUT(parentID,8) = segmentSetOUT(biggerChildID(1),8) + 1;
                                    
                                    tempRepFlag = 1;
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    [segmentSetOUT,notRepIntSizeErr] =...
                                        repairSizeFromIntensity(smallerChildID,segmentSetOUT,0,...
                                        notRepIntSizeErr,intensityInfo);
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                end

                            end                                
                             
                            clear tempChild1Chk tempChild2Chk
                            
                            notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'10',verbose);                             
                            
                        elseif (parentChk == 0)
                            %cases 2,5,8,11
                            tempRepFlag = 0;
                            if (isnan(segmentSetOUT(child1ID,3)) && isnan(segmentSetOUT(child2ID,3)))
                                %case 2
                                if ( (child1Size_rev == 1) || (child2Size_rev == 1) ) 
                                    segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                    if (checkSegmentSize(parentID,segmentSetOUT,-1,1,0) == 1)
                                        tempRepFlag = 1;
                                    else
                                        %Revert parent
                                        segmentSetOUT(parentID,8) = parentSize;
                                    end
                                end
                                
                                if (~tempRepFlag)                                    
                                    segmentSetOUT = assignSizeFromProb([parentID;NaN],parentSize,...
                                        [child1ID;child2ID],segmentSetIN,segmentSetOUT,intensityInfo);   
                                    
                                    tempRepFlag = 1;
                                end

                            elseif (isnan(segmentSetOUT(child1ID,3)) && ~isnan(segmentSetOUT(child2ID,3)) ) 
                                %case 5 = case 4                                
                                tempChild2Chk = checkSegmentSize(child2ID,segmentSetOUT,1,1,0);
                                tempRepFlag = 0;
                                if (tempChild2Chk == 1)
                                    %child ok.
                                    if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                        segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                        tempRepFlag = 1;
                                    elseif ((child1Size_rev > 1) && (child2Size_rev > 1)) 
                                        segmentSetOUT(child1ID,8) = 1;
                                        segmentSetOUT(parentID,8) = child2Size_rev + 1;
                                        tempRepFlag = 1;
                                    end
                                end
                                
                            elseif (~isnan(segmentSetOUT(child1ID,3)) && isnan(segmentSetOUT(child2ID,3)) ) 
                                %case 8 = case 7
                                tempChild1Chk = checkSegmentSize(child1ID,segmentSetOUT,1,1,0);
                                tempRepFlag = 0;
                                if (tempChild1Chk == 1)
                                    %child ok.
                                    if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                        segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                        tempRepFlag = 1;
                                    elseif ((child1Size_rev > 1) && (child2Size_rev > 1)) 
                                        segmentSetOUT(child2ID,8) = 1;
                                        segmentSetOUT(parentID,8) = child1Size_rev + 1;
                                        tempRepFlag = 1;
                                    end
                                end
                            elseif (~isnan(segmentSetOUT(child1ID,3)) && ~isnan(segmentSetOUT(child2ID,3)) ) 
                                %case 11 = case 10
                                tempChild1Chk = checkSegmentSize(child1ID,segmentSetOUT,1,1,0);
                                tempChild2Chk = checkSegmentSize(child2ID,segmentSetOUT,1,1,0);
                                tempRepFlag = 0;    
                                
                                if ( (tempChild1Chk == 1) && (tempChild2Chk == 1) )
                                    if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                        segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                        tempRepFlag = 1;                   
                                    else
                                        %Both > 1. The lower one must be set to
                                        %1 and then repaired.
                                        childIDs = [child1ID;child2ID];
                                        %childSizes = [child1Size_rev;child2Size_rev];
                                        childInts = [segmentSetIN(child1ID,7);segmentSetIN(child2ID,7)];
                                        %tempBigger = find((childInts == max(childInts)));
                                        smallerChildID = childIDs((childInts == min(childInts)));
                                        biggerChildID = childIDs(childIDs ~= smallerChildID);
                                        segmentSetOUT(smallerChildID(1),8) = 1;
                                        segmentSetOUT(parentID,8) = segmentSetOUT(biggerChildID(1),8) + 1;

                                        tempRepFlag = 1;
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        [segmentSetOUT,notRepIntSizeErr] =...
                                            repairSizeFromIntensity(smallerChildID,segmentSetOUT,0,...
                                            notRepIntSizeErr,intensityInfo);
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    end

                                end 
                            end

                            clear tempIntVals tempChild1Chk tempChild2Chk

                            notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'2,5,8,11',verbose);                                                         
                                
                        elseif (parentChk == 1)
                            %Both child segments continue
                            if (~isnan(segmentSetOUT(child1ID,3)) && ~isnan(segmentSetOUT(child2ID,3)) )
                                %case 12
                                tempChild1Chk = checkSegmentSize(child1ID,segmentSetOUT,1,1,0);
                                tempChild2Chk = checkSegmentSize(child2ID,segmentSetOUT,1,1,0);

                                tempRepFlag = 0;
                                %Can't repair if both segs have errors OR
                                %both are correct. Note: checkSegmentSize
                                %above can return NaN.                                    
                                if ( (tempChild1Chk == 1) && (tempChild2Chk == 1) )
                                    if (parentSize == 1)
                                        %parentChk = 1 is incorrect
                                        if ((child1Size_rev == 1) || (child2Size_rev == 1))
                                            segmentSetOUT(parentID,8) = child1Size_rev + child2Size_rev;
                                            tempRepFlag = 1;
                                        else
                                            %One child is also incorrect.
                                            %all three get updated.
                                            segmentSetOUT = assignSizeFromProb([parentID;NaN],parentSize,...
                                                [child1ID;child2ID],segmentSetIN,segmentSetOUT,intensityInfo); 
                                            tempRepFlag = 1;
                                            %need to recall all three.
                                        end
                                    else
                                        %Here only need to check if both
                                        %child segs > 1.
                                        if (child1Size_rev > 1) && (child2Size_rev > 1)
                                            %One child is also incorrect.
                                            %all three get updated.
                                            childIDs = [child1ID;child2ID];
                                            %childSizes = [child1Size_rev;child2Size_rev];
                                            childInts = [segmentSetIN(child1ID,7);segmentSetIN(child2ID,7)];
                                            %tempBigger = find((childInts == max(childInts)));
                                            smallerChildID = childIDs((childInts == min(childInts)));
                                            biggerChildID = childIDs(childIDs ~= smallerChildID);
                                            
                                            segmentSetOUT(smallerChildID(1),8) = 1; 
                                            
                                            %Smaller child set. Now, either
                                            %parent or bigger child or both
                                            %also need updating.
                                            %{
                                            segIntensity = [segmentSetIN(parentID,7);segmentSetIN(biggerChildID(1),7)];
                                            newParentSize = getSizeFromProb(parentSize,segIntensity,intensityInfo);
                                            segmentSetOUT(parentID,8) = newParentSize;
                                            segmentSetOUT(biggerChildID(1),8) = newParentSize - 1;
                                            %}
                                            if (parentSize == segmentSetOUT(biggerChildID(1),8) + 1)
                                                tempRepFlag = 1;
                                            else
                                                %Revert value to prevent
                                                %function call below.
                                                segmentSetOUT(smallerChildID(1),8) = ...
                                                    segmentSetIN(smallerChildID(1),8);
                                            end
                                            %need to call for child1 and 2.
                                        end
                                    end
                                                    
                                end                                           

                                notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                    [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'12',verbose);      
                            
                                if (segmentSetOUT(parentID,8) ~= parentSize)
                                    %Recall for parent with its parent
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    fprintf('\nRecall sp %d ',parentID);
                                    if (segmentSetOUT(parentID,10) == 1)
                                        segmentSetOUT(parentID,8) = parentSize;
                                        %tempRepFlag = 0;
                                    else
                                        segmentSetOUT(parentID,10) = 1;
                                        [segmentSetOUT,notRepIntSizeErr] =...
                                            repairSizeFromIntensity(segmentSetOUT(parentID,1),segmentSetOUT,0,...
                                            notRepIntSizeErr,intensityInfo);
                                        segmentSetOUT(parentID,10) = NaN;
                                    end
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                end

                                if ((tempChild1Chk == 1) &&...
                                        segmentSetOUT(child1ID,8) ~= child1Size_rev)
                                    %Recall with child1
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    fprintf('\nRecall c1 %d, parent %d ',child1ID,parentID);
                                    [segmentSetOUT,notRepIntSizeErr] =...
                                        repairSizeFromIntensity(child1ID,segmentSetOUT,0,...
                                        notRepIntSizeErr,intensityInfo);
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                end   
                                if ((tempChild2Chk == 1) &&...
                                        segmentSetOUT(child2ID,8) ~= child2Size_rev)
                                    %Recall with child2
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    fprintf('\nRecall c2 %d, parent %d',child2ID,parentID);
                                    [segmentSetOUT,notRepIntSizeErr] =...
                                        repairSizeFromIntensity(child2ID,segmentSetOUT,0,...
                                        notRepIntSizeErr,intensityInfo);
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                end                                  
                                

                                
                                clear tempChild1Chk tempChild2Chk tempRepFlag
                                
                            elseif (isnan(segmentSetOUT(child1ID,3)) && isnan(segmentSetOUT(child2ID,3)) )
                                %case 3
                                tempRepFlag = 0;
                                if (parentSize == 1)
                                    %parentChk = 1 is incorrect
                                    segmentSetOUT = assignSizeFromProb([parentID;NaN],parentSize,...
                                        [child1ID;child2ID],segmentSetIN,segmentSetOUT,intensityInfo); 
                                    
                                    tempRepFlag = 1;
                                elseif ( parentSize == 2 )
                                    %Neither child continues and parent's
                                    %size is 2.
                                    segmentSetOUT(child1ID,8) = 1;
                                    segmentSetOUT(child2ID,8) = 1;
                                    tempRepFlag = 1;
                                elseif (parentSize > 2)
                                    %If parent's size > 2, use existing
                                    %size values of the child segments to
                                    %repair. Not repaired if both equal.
                                    %Note: using intensity amplitude
                                    %because these values may differ for
                                    %two segments, but rounded they can be
                                    %equal.
                                    if (segmentSetOUT(child1ID,7) > segmentSetOUT(child2ID,7))
                                        segmentSetOUT(child1ID,8) = parentSize - 1;
                                        segmentSetOUT(child2ID,8) = 1;
                                        tempRepFlag = 1;
                                    elseif (segmentSetOUT(child1ID,7) < segmentSetOUT(child2ID,7))
                                        segmentSetOUT(child2ID,8) = parentSize - 1;
                                        segmentSetOUT(child1ID,8) = 1;
                                        tempRepFlag = 1;
                                    end
                                end
                                
                            notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'3',verbose);                                      
                                
                            elseif (isnan(segmentSetOUT(child1ID,3)) && ~isnan(segmentSetOUT(child2ID,3)) )
                                %case 6 ~ case 12 c
                                %First child ends, second child continues
                                tempChild2Chk = checkSegmentSize(child2ID,segmentSetOUT,1,1,0);
                                tempRepFlag = 0;                                
                                if (tempChild2Chk == 1)
                                    %child2 ok. 
                                    if (child2Size_rev > 1)
                                        segmentSetOUT(child1ID,8) = 1;
                                        if (parentSize == 1)
                                            %parentChk = 1 is incorrect
                                            segmentSetOUT(parentID,8) = child2Size_rev + 1;
                                            tempRepFlag = 1;
                                        else
                                            if ( parentSize == child2Size_rev + 1 ) 
                                                tempRepFlag = 1;
                                            else
                                                %parent has a mismatch with child
                                                %segments - parent and child 2
                                                %are ok, with child 1 set to 1,
                                                %but the sum is not working. so
                                                %child 1 and 2 are updated,
                                                %parent kept as is.
                                                childIDs = [child1ID;child2ID];
                                                childInts = [segmentSetIN(child1ID,7);segmentSetIN(child2ID,7)];
                                                smallerChildID = childIDs((childInts == min(childInts)));
                                                biggerChildID = childIDs(childIDs ~= smallerChildID);
                                                segmentSetOUT(smallerChildID(1),8) = 1; 
                                                segmentSetOUT(biggerChildID(1),8) = ...
                                                    parentSize - 1;
                                                %Need to recall if child2 size
                                                %has been changed
                                                if (segmentSetOUT(child2ID,8) ~= child2Size_rev)
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    [segmentSetOUT,notRepIntSizeErr] =...
                                                        repairSizeFromIntensity(child2ID,segmentSetOUT,0,...
                                                        notRepIntSizeErr,intensityInfo);
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                end

                                                tempRepFlag = 1;
                                            end
                                        end
                                    elseif (child2Size_rev == 1)
                                        if (parentSize == 1)
                                            %parentChk = 1 is incorrect
                                            segIntensity = [segmentSetIN(parentID,7);segmentSetIN(child1ID,7)];
                                            newParentSize = getSizeFromProb(parentSize,segIntensity,intensityInfo);
                                            segmentSetOUT(parentID,8) = newParentSize;
                                            segmentSetOUT(child1ID,8) = newParentSize - 1;
                                            tempRepFlag = 1;
                                        else
                                            segmentSetOUT(child1ID,8) = parentSize - 1;
                                            tempRepFlag = 1;
                                        end
                                    end
                                        
                                end    
    
                                %021414
                                if (segmentSetOUT(parentID,8) ~= parentSize)
                                    if (isnan(segmentSetOUT(parentID,10)))
                                        segmentSetOUT(parentID,10) = 1;
                                        %Recall for parent with its parent
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        [segmentSetOUT,notRepIntSizeErr] =...
                                            repairSizeFromIntensity(segmentSetOUT(parentID,1),segmentSetOUT,0,...
                                            notRepIntSizeErr,intensityInfo);
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    else
                                        segmentSetOUT(parentID,10) = NaN;
                                        tempRepFlag = 0;
                                    end
                                end                                              
                                
                                clear tempChild2Chk
                                
                                notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                    [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'6',verbose);                                      
                                
                            elseif (~isnan(segmentSetOUT(child1ID,3)) && isnan(segmentSetOUT(child2ID,3)) )
                                %case 9 ~ case 12d
                                tempChild1Chk = checkSegmentSize(child1ID,segmentSetOUT,1,1,0);
                                tempRepFlag = 0;
                                if (tempChild1Chk == 1)
                                    %child1 ok. 
                                    if (child1Size_rev > 1)
                                        segmentSetOUT(child2ID,8) = 1;
                                        if (parentSize == 1)
                                            %parentChk = 1 is incorrect
                                            segmentSetOUT(parentID,8) = child1Size_rev + 1;
                                            tempRepFlag = 1;
                                        else
                                            if ( parentSize == child1Size_rev + 1) 
                                                tempRepFlag = 1;
                                            else                                                
                                                %there is a mismatch with child
                                                %segments
                                                childIDs = [child1ID;child2ID];
                                                childInts = [segmentSetIN(child1ID,7);segmentSetIN(child2ID,7)];
                                                smallerChildID = childIDs((childInts == min(childInts)));
                                                biggerChildID = childIDs(childIDs ~= smallerChildID);
                                                segmentSetOUT(smallerChildID(1),8) = 1; 
                                                segmentSetOUT(biggerChildID(1),8) = ...
                                                    parentSize - 1;
                                                %Need to recall if child2 size
                                                %has been changed
                                                if (segmentSetOUT(child1ID,8) ~= child1Size_rev)
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    [segmentSetOUT,notRepIntSizeErr] =...
                                                        repairSizeFromIntensity(child1ID,segmentSetOUT,0,...
                                                        notRepIntSizeErr,intensityInfo);
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                end
                                                tempRepFlag = 1;
                                            end
                                        end
                                    elseif (child1Size_rev == 1)
                                        if (parentSize == 1)
                                            %parentChk = 1 is incorrect
                                            segIntensity = [segmentSetIN(parentID,7);segmentSetIN(child1ID,7)];
                                            newParentSize = getSizeFromProb(parentSize,segIntensity,intensityInfo);
                                            segmentSetOUT(parentID,8) = newParentSize;
                                            segmentSetOUT(child2ID,8) = newParentSize - 1;
                                            tempRepFlag = 1;
                                        else
                                            segmentSetOUT(child2ID,8) = parentSize - 1;
                                            tempRepFlag = 1;
                                        end
                                    end
                                end %tempChild1Chk

                                %022414
                                if (segmentSetOUT(parentID,8) ~= parentSize)
                                    if (isnan(segmentSetOUT(parentID,10)))
                                        segmentSetOUT(parentID,10) = 1;
                                        %Recall for parent with its parent
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        [segmentSetOUT,notRepIntSizeErr] =...
                                            repairSizeFromIntensity(segmentSetOUT(parentID,1),segmentSetOUT,0,...
                                            notRepIntSizeErr,intensityInfo);
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    else
                                        segmentSetOUT(parentID,10) = NaN;
                                        tempRepFlag = 0;
                                    end
                                end 

                                clear tempChild1Chk
                                
                                notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parentID;NaN],...
                                    [child1ID;child2ID],segmentSetOUT,notRepIntSizeErr,'9',verbose);                                    
                            end

                        end %parentChk

                    end %if ~swapFxd
                                
                end %if child 2 revised

            else
                fprintf('Repaired (s3): (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,parentSize,...
                    child1ID,child1Size_rev,child2ID,child2Size);                    
            end %if child 1 revised

        elseif (checkAll && (child1Size + child2Size == parentSize))
            %Continue on with child 1
            [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(segmentSetIN(segID,3),...
                segmentSetIN,1,notRepIntSizeErr,intensityInfo);
            %Continue on with child 2
            [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(segmentSetOUT(segID,4),...
                segmentSetIN,1,notRepIntSizeErr,intensityInfo);            
        else
            %No problem found.
            %{
            fprintf(' OK: (%d:%d) -> (%d:%d) + (%d:%d). ',parentID,parentSize,...
                child1ID,child1Size,child2ID,child2Size);                                
            %}
        end %if mismatch exists
        
    elseif (~isnan(segmentSetIN(segID,3)) && ~isnan(segmentSetIN(segID,6)) )
        %Merge event.
        parent1ID = segID;
        parent2ID = segmentSetIN(parent1ID,6);
        childID = segmentSetIN(parent1ID,3);
        parent1Size = segmentSetIN(parent1ID,8);
        parent2Size = segmentSetIN(parent2ID,8);
        childSize = segmentSetIN(childID,8);

        if ( (childSize ~= parent1Size + parent2Size) ||...
                (childSize == 0) || (parent1Size == 0) || (parent2Size == 0) ||...
                ((parent1Size > 1) && (parent2Size > 1)) )
            %{
            fprintf('\n(%d:%d) + (%d:%d) -> (%d:%d). ',segID,parent1Size,...
                parent2ID,parent2Size,childID,childSize);
            %}

            [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(childID,segmentSetIN,0,notRepIntSizeErr,intensityInfo);
            childSize_rev = segmentSetOUT(childID,8);                    
            %Will use segmentSetOUT instead of IN so that the child's value
            %is used throughout. The value returned above is either the
            %same as original or has been updated to fix subsequent events.
            %Thus, keep value.
            if ( (childSize_rev ~= parent1Size + parent2Size) || (childSize_rev == 0) || ...
                (parent1Size == 0) || (parent2Size == 0) ||...
                ((parent1Size > 1) && (parent2Size > 1)) )               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Try to repair via swap
                swapFxd = 0;
                if ( (childSize_rev > 0) && (parent1Size > 0) && (parent2Size > 0))
                    parent1SibSize = NaN;
                    parent2SibSize = NaN;
                    childPartSize = NaN;
                    %Get child's partner's size
                    if (~isnan(segmentSetOUT(childID,6)) && isnan(segmentSetOUT(segmentSetOUT(childID,6),1)) )
                        childPartSize = segmentSetOUT(segmentSetOUT(childID,6),8);
                    end
                    %Get parent1's sibling's size
                    if (~isnan(segmentSetOUT(parent1ID,5)) && isnan(segmentSetOUT(segmentSetOUT(parent1ID,5),3)) )
                        parent1SibSize = segmentSetOUT(segmentSetOUT(parent1ID,5),8);
                    end
                    %Get parent2's sibling's size
                    %NOTE: parent1 and parent2 maybe eachother's siblings - ADD
                    if (~isnan(segmentSetOUT(parent2ID,5)) && isnan(segmentSetOUT(segmentSetOUT(parent2ID,5),3)) )
                        parent2SibSize = segmentSetOUT(segmentSetOUT(parent2ID,5),8);
                    end
                    %All possible sums with swapped values
                    parentVec1 = [parent1Size; parent1Size; parent1SibSize; parent1SibSize];
                    parentVec2 = [parent2Size; parent2SibSize; parent2Size; parent2SibSize];
                    parentSumVec = parentVec1 + parentVec2;
                    correctSum = find(parentSumVec == childPartSize,1,'first');
                    %NOTE: can have more than one matching sum
                    if (~isempty(correctSum))
                        %A swap that works found. Assign values.
                        %Child takes sibling's value and vice versaw
                        segmentSetOUT(segmentSetOUT(childID,6),8) = childSize_rev;
                        segmentSetOUT(childID,8) = childPartSize;
                        %Parent 1
                        if (parentVec1(correctSum) == parent1SibSize)
                            segmentSetOUT(segmentSetOUT(parent1ID,5),8) = parent1Size;    
                            segmentSetOUT(parent1ID,8) = parentVec1(correctSum);
                        end
                        %Parent 2
                        if (parentVec2(correctSum) == parent2SibSize)
                            segmentSetOUT(segmentSetOUT(parent2ID,5),8) = parent2Size;    
                            segmentSetOUT(parent2ID,8) = parentVec2(correctSum);
                        end
                        swapFxd = 1;
                        fprintf('Repaired (swap ma): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                            parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));
                    else
                        correctSum = find(parentSumVec(2:4) == childSize_rev,1,'first');
                        if (correctSum)
                            %Only swapping parents here.
                            if (parentVec1(correctSum) == parent1SibSize)
                                segmentSetOUT(segmentSetOUT(parent1ID,5),8) = parent1Size;    
                                segmentSetOUT(parent1ID,8) = parentVec1(correctSum);
                            end
                            if (parentVec2(correctSum) == parent2SibSize)
                                segmentSetOUT(segmentSetOUT(parent2ID,5),8) = parent2Size;    
                                segmentSetOUT(parent2ID,8) = parentVec2(correctSum);
                            end
                            swapFxd = 1;    
                            fprintf('Repaired (swap mb): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                                parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));
                        end

                    end %Size check 1
                end %swapfx
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (~swapFxd)
                    %Attempt to repair via assignment
                    %First, check if parents' prior events are ok. 
                    parent1Chk = checkSegmentSize(parent1ID,segmentSetOUT,-1,1,0);
                    parent2Chk = checkSegmentSize(parent2ID,segmentSetOUT,-1,1,0);
                    
                    if (isnan(parent1Chk) && isnan(parent2Chk) && isnan(segmentSetIN(childID,3)))
                        %case 1
                        segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                            [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);                         
                        %{
                        outStr = sprintf('Not repaired (1): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                        parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));
                        errMagnitude = (segmentSetOUT(parent1ID,8) + ...
                                segmentSetOUT(parent2ID,8)) - segmentSetOUT(childID,8);                           
                        notRepIntSizeErr = setNotRepErrList(notRepIntSizeErr,[parent1ID;parent2ID],'1',...
                            outStr,'merge',errMagnitude);
                        fprintf(outStr);                        
                        %}
                    elseif (isnan(parent1Chk) && isnan(parent2Chk) && ~isnan(segmentSetIN(childID,3))) 
                        %case 10
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);
                        tempRepFlag = 0;
                        if (tempChildChk ~= 1)
                            fprintf('\nChild not ok found (10).');
                        else
                            %Child is ok.
                            if (childSize_rev == 2)
                                segmentSetOUT(parent1ID,8) = 1;
                                segmentSetOUT(parent2ID,8) = 1;
                                tempRepFlag = 1;
                            elseif (childSize_rev > 2)                                
                                if (segmentSetOUT(parent1ID,7) > segmentSetOUT(parent2ID,7))
                                    segmentSetOUT(parent1ID,8) = childSize_rev - 1;
                                    segmentSetOUT(parent2ID,8) = 1;
                                    tempRepFlag = 1;
                                elseif (segmentSetOUT(parent1ID,7) < segmentSetOUT(parent2ID,7))
                                    segmentSetOUT(parent2ID,8) = childSize_rev - 1;
                                    segmentSetOUT(parent1ID,8) = 1;
                                    tempRepFlag = 1;
                                end
                            elseif (childSize_rev == 1)
                                %ChildChk ok is incorrect.
                                segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                    [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);  

                                tempRepFlag = 1; 
                                %if childChk = 1 is incorrect
                                %if (segmentSetOUT(childID,8) ~= childSize_rev)
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    [segmentSetOUT,notRepIntSizeErr] =...
                                        repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                        notRepIntSizeErr,intensityInfo);
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                %end                                   
                            end

                        end
                     
                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'10',verbose);                             
                            
                    elseif ((parent1Chk == 0) && (parent2Chk == 0))
                        tempRepFlag = 0;
                        if (isnan(segmentSetIN(childID,3)))
                            %case 5
                            segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);                             
                           
                            tempRepFlag = 1;
                        else
                            %case 14
                            if (checkSegmentSize(childID,segmentSetOUT,1,1,0) == 1)
                                %case 14 with child ok.
                                %try updating either parent 
                                if (childSize_rev > 1)
                                    %Set smaller one to 1. Set other to
                                    %parent - 1. It may or may not give
                                    %parent1Chk = 1 and/or parent2Chk = 1.                                    
                                    parentIDs = [parent1ID;parent2ID];                            
                                    parentInts = [segmentSetIN(parent1ID,7);segmentSetIN(parent2ID,7)];                            
                                    biggerParentID = parentIDs((parentInts == max(parentInts)));
                                    smallerParentID = parentIDs(parentIDs ~= biggerParentID(1));                                    
                                    
                                    segmentSetOUT(smallerParentID(1),8) = 1; 
                                    segmentSetOUT(biggerParentID(1),8) = childSize_rev - 1;

                                    tempRepFlag = 1;                                                                                               
                                else
                                    %child not ok.
                                    segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                        [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);  
                                    
                                    tempRepFlag = 1;   

                                    %if childChk = 1 is incorrect
                                    %if (segmentSetOUT(childID,8) ~= childSize_rev)
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        [segmentSetOUT,notRepIntSizeErr] =...
                                            repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                            notRepIntSizeErr,intensityInfo);
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                    %end                                     
                                end
                                
                            end


                        end                   
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'5,14',verbose);                             
                        
                    elseif ((parent1Chk == 1) && (parent2Chk == 1))
                        if (isnan(segmentSetIN(childID,3)))
                            %case 9
                            if ( ((parent1Size == 1) || (parent2Size == 1)) )
                                segmentSetOUT(childID,8) = parent1Size + parent2Size;
                                fprintf('Repaired (9): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                                    parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));                                
                            else
                                %Both parent sizes > 1.
                                %1. The lower one must be set to
                                %1 and then repaired.  (split 10a)
                                %Need new value for one parent and new
                                %value for child.
                                parentIDs = [parent1ID;parent2ID];                            
                                parentInts = [segmentSetIN(parent1ID,7);segmentSetIN(parent2ID,7)];                            
                                biggerParentID = parentIDs((parentInts == max(parentInts)));
                                smallerParentID = parentIDs(parentIDs ~= biggerParentID(1));
                                segmentSetOUT(smallerParentID(1),8) = 1;  
                                
                                segmentSetOUT(childID,8) = ...
                                    segmentSetOUT(biggerParentID(1),8) + 1;
                                    
                                %Since parentChk = 1 was incorrect
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if (isnan(segmentSetOUT(smallerParentID(1),10)))
                                    segmentSetOUT(smallerParentID(1),10) = 1; 
                                    %Recall for parent with its parent
                                    [segmentSetOUT,notRepIntSizeErr] =...
                                        repairSizeFromIntensity(segmentSetOUT(smallerParentID(1),1),segmentSetOUT,0,...
                                        notRepIntSizeErr,intensityInfo);
                                else
                                    segmentSetOUT(smallerParentID(1),10) = NaN;
                                    %tempRepFlag = 0;
                                end                                    

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                   
                                

                                %{
                                outStr = sprintf('Not repaired (9): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                                    parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));
                                errMagnitude = (segmentSetOUT(parent1ID,8) + ...
                                        segmentSetOUT(parent2ID,8)) - segmentSetOUT(childID,8);                           
                                notRepIntSizeErr = setNotRepErrList(notRepIntSizeErr,[parent1ID;parent2ID],'9',...
                                    outStr,'merge',errMagnitude);
                                fprintf(outStr);
                                %}
                            end                                
                        else
                            %case 18
                            tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);
                            tempRepFlag = 0;
                            if (tempChildChk ~= 1)
                               fprintf('\nChild not ok found (18).');
                            else
                                %child ok.
                                %if childSize_rev > 1 and (p1Size = 1 or
                                %p2Size = 1), not repairable. Else...
                                if ((parent1Size == 1) || (parent2Size == 1))
                                    if (childSize_rev == 1)
                                        %childChk = 1 is incorrect
                                        segmentSetOUT(childID,8) =...
                                            segmentSetOUT(parent1ID,8) + segmentSetOUT(parent2ID,8);
                                        tempRepFlag = 1;
                                    end
                                    %if childSize_rev > 1, no repair.
                                elseif ((parent1Size > 1) && (parent2Size > 1))
                                    %Set smaller parent to 1 and recall.
                                    parentIDs = [parent1ID;parent2ID];                            
                                    parentInts = [segmentSetIN(parent1ID,7);segmentSetIN(parent2ID,7)];                            
                                    biggerParentID = parentIDs((parentInts == max(parentInts)));
                                    smallerParentID = parentIDs(parentIDs ~= biggerParentID(1));
                                    segmentSetOUT(smallerParentID(1),8) = 1;  
                                    %{
                                    segmentSetOUT(childID,8) = ...
                                        segmentSetOUT(biggerParentID(1),8) + 1;
                                    %}
                                    if ((childSize_rev > 1) &&...
                                            (childSize_rev ==...
                                            segmentSetOUT(biggerParentID(1),8) + 1))
                                        tempRepFlag = 1;
                                    elseif (childSize_rev == 1)
                                        %childChk = 1 is incorrect
                                        segmentSetOUT(childID,8) =...
                                            segmentSetOUT(parent1ID,8) + segmentSetOUT(parent2ID,8);
                                        tempRepFlag = 1;
                                    else
                                        %Revert to prevent function call.
                                        segmentSetOUT(smallerParentID(1),8) =...
                                            segmentSetIN(smallerParentID(1),8);
                                    end
                                        
                                    %If parentChk = 1 was incorrect
                                   if (segmentSetOUT(smallerParentID(1),8) ~=...
                                            segmentSetIN(smallerParentID(1),8))
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        if (isnan(segmentSetOUT(smallerParentID(1),10)))
                                            segmentSetOUT(smallerParentID(1),10) = 1;
                                            %Recall for parent with its parent
                                            [segmentSetOUT,notRepIntSizeErr] =...
                                                repairSizeFromIntensity(segmentSetOUT(smallerParentID(1),1),segmentSetOUT,0,...
                                                notRepIntSizeErr,intensityInfo);
                                        else
                                            segmentSetOUT(smallerParentID(1),10) = NaN;
                                            tempRepFlag = 0;
                                        end
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                       
                                   end
                                end

                                %If childChk = 1 was incorrect or child
                                %changed - only in this block
                                if (segmentSetOUT(childID,8) ~= childSize_rev)
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    [segmentSetOUT,notRepIntSizeErr] =...
                                        repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                        notRepIntSizeErr,intensityInfo);
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                     
                                end
                            end

                            clear tempChildChk
                            
                            notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                                childID,segmentSetOUT,notRepIntSizeErr,'18',verbose);                              
                            
                        end %if isnan child
                    elseif ( ((isnan(parent1Chk) && (parent2Chk == 1)) || ...
                                ((parent1Chk == 1) && isnan(parent2Chk)) ) && ...
                                    isnan(segmentSetIN(childID,3)) )                                
                                %cases 3 & 7  
                                if (isnan(parent1Chk))
                                    %case 3
                                    if (parent2Size == 1)
                                        segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent1ID,7)];
                                        newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                        segmentSetOUT(childID,8) = newChildSize;
                                        segmentSetOUT(parent1ID,8) = newChildSize - 1;  
                                    elseif (parent2Size > 1)
                                        segmentSetOUT(parent1ID,8) = 1;
                                        segmentSetOUT(childID,8) = parent2Size + 1;                                        
                                    end
                                elseif (isnan(parent2Chk))
                                    %case 7
                                    if (parent1Size == 1)
                                        segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent2ID,7)];
                                        newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                        segmentSetOUT(childID,8) = newChildSize;
                                        segmentSetOUT(parent2ID,8) = newChildSize - 1;     
                                    elseif (parent1Size > 1)
                                        segmentSetOUT(parent2ID,8) = 1;
                                        segmentSetOUT(childID,8) = parent1Size + 1;                                        
                                    end
                                end
                            %{                            
                            outStr = sprintf('Not repaired (3,7): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                                parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));  
                            errMagnitude = (segmentSetOUT(parent1ID,8) + ...
                                    segmentSetOUT(parent2ID,8)) - segmentSetOUT(childID,8);                           
                            notRepIntSizeErr = setNotRepErrList(notRepIntSizeErr,[parent1ID;parent2ID],'12',...
                                outStr,'merge',errMagnitude);  
                            fprintf(outStr);   
                            %}
                    elseif ( isnan(parent1Chk) && (parent2Chk == 1) )
                        %case 12 ~ case 15
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);
                        tempRepFlag = 0;
                        if ( tempChildChk ~= 1 )
                            fprintf('\nChild not ok found (12).');
                        else
                            %child ok.
                            if (parent2Size > 1)
                                segmentSetOUT(parent1ID,8) = 1;
                                if ( (childSize_rev == 1) ||...
                                        (parent2Size + 1 ~= childSize_rev))
                                    %there is a mismatch with child
                                    segmentSetOUT(childID,8) = parent2Size + 1;                                    
                                    %need to ReCall with child
                                end
                                tempRepFlag = 1;
                            elseif (parent2Size == 1)
                                if (childSize_rev > 1)
                                    segmentSetOUT(parent1ID,8) = childSize_rev - 1;                                    
                                elseif (childSize_rev == 1)
                                    %childChk = 1 is incorrect
                                    segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent1ID,7)];
                                    newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                    segmentSetOUT(childID,8) = newChildSize;
                                    segmentSetOUT(parent1ID,8) = newChildSize - 1;
                                end
                                tempRepFlag = 1;                                    
                            end       
                  
                            if (segmentSetOUT(childID,8) ~= childSize_rev)
                                %Since childChk = 1 was incorrect
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                [segmentSetOUT,notRepIntSizeErr] =...
                                    repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                    notRepIntSizeErr,intensityInfo);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                            end 
                                
                        end
  
                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'12',verbose);                             
                        
                    elseif ( (parent1Chk == 1) && isnan(parent2Chk) )
                        %case 16 ~ case 17
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);
                        tempRepFlag = 0;
                        if ( tempChildChk ~= 1 )
                            fprintf('\nChild not ok found (16).');
                        else
                            %child ok.
                           if (parent1Size > 1)
                                segmentSetOUT(parent2ID,8) = 1;
                                if ( (childSize_rev == 1) ||...
                                        (parent1Size + 1 ~= childSize_rev))
                                    %there is a mismatch with child
                                    segmentSetOUT(childID,8) = parent1Size + 1;                                    
                                    %need to ReCall with child
                                end
                                tempRepFlag = 1;
                            elseif (parent1Size == 1)
                                if (childSize_rev > 1)
                                    segmentSetOUT(parent2ID,8) = childSize_rev - 1;                                    
                                elseif (childSize_rev == 1)
                                    %childChk = 1 is incorrect
                                    segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent2ID,7)];
                                    newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                    segmentSetOUT(childID,8) = newChildSize;
                                    segmentSetOUT(parent2ID,8) = newChildSize - 1;
                                end
                                tempRepFlag = 1;                                    
                            end       

                            %Since childChk = 1 was incorrect
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            [segmentSetOUT,notRepIntSizeErr] =...
                                repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                notRepIntSizeErr,intensityInfo);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                        end
 
                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'16',verbose);                                                     
                        
                    elseif ( (parent1Chk ~= 0) && (parent2Chk == 0) && ...
                            isnan(segmentSetIN(childID,3)))
                        %case 2, 8
                        tempRepFlag = 0;
                        if (isnan(parent1Chk))
                            %case 2 = case 4
                            if (childSize_rev > 1)                                
                                if ( (childSize_rev > parent1Size) &&...
                                        (parent1Size > 1) || (childSize_rev - parent1Size > 1))                                
                                    segmentSetOUT(parent2ID,8) = childSize_rev - parent1Size;
                                    %Recheck
                                    if(checkSegmentSize(parent2ID,segmentSetOUT,-1,1,0) == 1)
                                        tempRepFlag = 1;                                
                                    end
                                end
                                if (~tempRepFlag)
                                    %Update failed. Revert value.
                                    segmentSetOUT(parent2ID,8) = parent2Size;

                                    segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                        [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);                                      

                                    tempRepFlag = 1;                                                              
                                end                                
                            else
                                segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                    [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);  

                                tempRepFlag = 1;   
                            end
                        elseif (parent1Chk == 1)
                            %case 8 = case 6
                            if (parent1Size == 1)
                                %Need to update parent2 and child
                                segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent2ID,7)];
                                newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                segmentSetOUT(childID,8) = newChildSize;
                                segmentSetOUT(parent2ID,8) = newChildSize - 1;
                                tempRepFlag = 1;                                 
                            elseif (parent1Size > 1)
                                segmentSetOUT(parent2ID,8) = 1;
                                segmentSetOUT(childID,8) = parent1Size + 1;
                                tempRepFlag = 1;
                            end                           
                        end
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'2,8',verbose);                             
                        
                    elseif ( isnan(parent1Chk) && (parent2Chk == 0) && ...
                            ~isnan(segmentSetIN(childID,3)))
                        %case 11 = case 13
                        tempRepFlag = 0;
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);
                        if (tempChildChk == 1) 
                            %child ok
                            parentIDs = [parent1ID;parent2ID];                            
                            parentInts = [segmentSetIN(parent1ID,7);segmentSetIN(parent2ID,7)];                            
                            biggerParentID = parentIDs((parentInts == max(parentInts)));
                            smallerParentID = parentIDs(parentIDs ~= biggerParentID(1));
                            %Set smaller parent to 1.
                            segmentSetOUT(smallerParentID(1),8) = 1;  
                            
                            if (childSize_rev > 1)
                                segmentSetOUT(biggerParentID(1),8) =...
                                    childSize_rev - 1;
                                tempRepFlag = 1;
                            elseif (childSize_rev == 1)                                    
                                %childChk ok is incorrect. need to also
                                %update child.
                                segIntensity = [segmentSetIN(childID,7);segmentSetIN(biggerParentID(1),7)];
                                newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);
                                
                                segmentSetOUT(childID,8) = newChildSize;
                                segmentSetOUT(biggerParentID(1),8) =...
                                    newChildSize - 1;                                
                                %Since childChk = 1 was incorrect
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                [segmentSetOUT,notRepIntSizeErr] =...
                                    repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                    notRepIntSizeErr,intensityInfo);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                tempRepFlag = 1;
                            end
                        else
                            fprintf('\nChild not ok found (11).');
                        end
                                             
                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'11',verbose);                                                                             
                        
                    elseif ( (parent1Chk == 1) && (parent2Chk == 0) && ...
                            ~isnan(segmentSetIN(childID,3)))
                        %case 17 ~ case 16
                        tempRepFlag = 0;
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);                                                                  
                        if ( tempChildChk ~= 1 )
                            fprintf('\nChild not ok found (17).');
                        else
                            %child ok.
                            if (parent1Size > 1)
                                segmentSetOUT(parent2ID,8) = 1;
                                if ( (childSize_rev == 1) ||...
                                        (parent1Size + 1 ~= childSize_rev))
                                    %there is a mismatch with child
                                    segmentSetOUT(childID,8) = parent1Size + 1;                                    
                                    %need to ReCall with child
                                end
                                tempRepFlag = 1;
                            elseif (parent1Size == 1)
                                if (childSize_rev > 1)
                                    segmentSetOUT(parent2ID,8) = childSize_rev - 1;                                    
                                elseif (childSize_rev == 1)
                                    %childChk = 1 is incorrect
                                    segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent2ID,7)];
                                    newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                    segmentSetOUT(childID,8) = newChildSize;
                                    segmentSetOUT(parent2ID,8) = newChildSize - 1;
                                end
                                tempRepFlag = 1;                                    
                            end     
                            if (segmentSetOUT(childID,8) ~= childSize_rev)
                                %Since childChk = 1 was incorrect
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                [segmentSetOUT,notRepIntSizeErr] =...
                                    repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                    notRepIntSizeErr,intensityInfo);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                            end 
                            
                        end  
                        
                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'17',verbose);                             
                        
                    elseif ( (parent1Chk == 0) && (parent2Chk ~= 0) && ...
                            isnan(segmentSetIN(childID,3)))
                        tempRepFlag = 0;
                        if (isnan(parent2Chk))
                            %case 4 = case 2
                            if (childSize_rev > 1)                               
                                if ( (childSize_rev > parent2Size) &&...
                                        (parent2Size > 1) || (childSize_rev - parent2Size > 1))                                
                                    segmentSetOUT(parent1ID,8) = childSize_rev - parent2Size;
                                    %Recheck
                                    if(checkSegmentSize(parent1ID,segmentSetOUT,-1,1,0) == 1)
                                        tempRepFlag = 1;                                
                                    end
                                end
                                if (~tempRepFlag)
                                    %Update failed. Revert value.
                                    segmentSetOUT(parent1ID,8) = parent1Size;
 
                                    segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                        [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);                                      

                                    tempRepFlag = 1;                                                              
                                end                                
                            else
                                
                                segmentSetOUT = assignSizeFromProb([parent1ID;parent2ID],childSize_rev,...
                                    [childID;NaN],segmentSetIN,segmentSetOUT,intensityInfo);  

                                tempRepFlag = 1;   
                            end
                        elseif (parent2Chk == 1)
                            %case 6 = case 8
                            if (parent2Size == 1)
                                %Need to update parent1 and child
                                segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent1ID,7)];
                                newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                segmentSetOUT(childID,8) = newChildSize;
                                segmentSetOUT(parent1ID,8) = newChildSize - 1;
                                tempRepFlag = 1;                                 
                            elseif (parent2Size > 1)
                                segmentSetOUT(parent1ID,8) = 1;
                                segmentSetOUT(childID,8) = parent2Size + 1;
                                tempRepFlag = 1;
                            end                           
                        end                       

                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'4,6',verbose);                                                     

                    elseif ( (parent1Chk == 0) && isnan(parent2Chk) && ...
                            ~isnan(segmentSetIN(childID,3)))                        
                        %case 13 = case 11
                        tempRepFlag = 0;
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);
                        if (tempChildChk == 1) 
                            %child ok
                            parentIDs = [parent1ID;parent2ID];                            
                            parentInts = [segmentSetIN(parent1ID,7);segmentSetIN(parent2ID,7)];                            
                            biggerParentID = parentIDs((parentInts == max(parentInts)));
                            smallerParentID = parentIDs(parentIDs ~= biggerParentID(1));
                            %Set smaller parent to 1.
                            segmentSetOUT(smallerParentID(1),8) = 1;
                            if (childSize_rev > 1)
                                segmentSetOUT(biggerParentID(1),8) =...
                                    childSize_rev - 1;
                                tempRepFlag = 1;                                
                            elseif (childSize_rev == 1)                                    
                                %childChk ok is incorrect. need to also
                                %update child.
                                segIntensity = [segmentSetIN(childID,7);segmentSetIN(biggerParentID(1),7)];
                                newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);
                                
                                segmentSetOUT(childID,8) = newChildSize;
                                segmentSetOUT(biggerParentID(1),8) =...
                                    newChildSize - 1;                                
                                %Since childChk = 1 was incorrect
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                [segmentSetOUT,notRepIntSizeErr] =...
                                    repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                    notRepIntSizeErr,intensityInfo);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                                tempRepFlag = 1;
                            end
                        else
                            fprintf('\nChild not ok found (13).');
                        end                        

                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'13',verbose);                                                     
                        
                    elseif ( (parent1Chk == 0) && (parent2Chk == 1) && ...
                            ~isnan(segmentSetIN(childID,3)))                        
                        %case 15 = case 17
                        tempRepFlag = 0;
                        tempChildChk = checkSegmentSize(childID,segmentSetOUT,1,1,0);                                                                  
                        if ( tempChildChk ~= 1 )
                            %child not ok.                            
                            fprintf('\nChild not ok found (15).');                            
                        else
                            %child ok.
                            if (parent2Size > 1)
                                segmentSetOUT(parent1ID,8) = 1;
                                if ( (childSize_rev == 1) ||...
                                        (parent2Size + 1 ~= childSize_rev))
                                    %there is a mismatch with child
                                    segmentSetOUT(childID,8) = parent2Size + 1;                                    
                                    %need to ReCall with child
                                end
                                tempRepFlag = 1;
                            elseif (parent2Size == 1)
                                if (childSize_rev > 1)
                                    segmentSetOUT(parent1ID,8) = childSize_rev - 1;                                    
                                elseif (childSize_rev == 1)
                                    %childChk = 1 is incorrect
                                    segIntensity = [segmentSetIN(childID,7);segmentSetIN(parent1ID,7)];
                                    newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

                                    segmentSetOUT(childID,8) = newChildSize;
                                    segmentSetOUT(parent1ID,8) = newChildSize - 1;
                                end
                                tempRepFlag = 1;                                    
                            end     
                            if (segmentSetOUT(childID,8) ~= childSize_rev)
                                %Since childChk = 1 was incorrect
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                [segmentSetOUT,notRepIntSizeErr] =...
                                    repairSizeFromIntensity(childID,segmentSetOUT,0,...
                                    notRepIntSizeErr,intensityInfo);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                            end   
                        end  

                        clear tempChildChk
                        
                        notRepIntSizeErr = perCaseRepairResult(tempRepFlag,[parent1ID;parent2ID],...
                            childID,segmentSetOUT,notRepIntSizeErr,'15',verbose);                             
                        
                    end

                end %if ~swapfxd
            
                
            else %On merge, if revised child value did not work.
                %Revised value worked.
                fprintf('Repaired (m1): (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                    parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8));                
            end            
        elseif (checkAll && (childSize == parent1Size + parent2Size))
            %Continue on with child 1
            [segmentSetOUT,notRepIntSizeErr] = repairSizeFromIntensity(segmentSetIN(segID,3),...
                segmentSetIN,1,notRepIntSizeErr,intensityInfo);            
        else %if mismatch exists
            %{
            fprintf(' OK: (%d:%d) + (%d:%d) -> (%d:%d). ',parent1ID,segmentSetOUT(parent1ID,8),...
                    parent2ID,segmentSetOUT(parent2ID,8),childID,segmentSetOUT(childID,8)); 
            %}
        end
        
    end %main if
    
    
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%BEGIN SUBFUNCTIONS

function segmentSetOUT = assignSizeFromProb(parentID,segSize,childID,...
    segmentSetIN,segmentSetOUT,intensityInfo)
%ASSIGNSIZEFROMPROB uses segments' intensity amplitude, the currently
%determined sizes and the intensity quantum to probablistically
%assign the most likely sizes. This approach is used for cases where
%determining sizes based on events is not possible.
%   
%   INPUT:
%           parentID:           segment ID of parents (1 or 2 element)
%           segSize:            size of parent or child determined but to
%                               be adjusted here possible
%           segmetSetIN:        original table of segments and events
%                               constructed by processSegmentEvents
%           segmetSetOUT:       current table of segments and events
%                               constructed by processSegmentEvents
%           intensityInfo:      intensity quantum as defined in main sim.
%
%   OUTPUT:
%           segmetSetOUT:       updated table of segments and events
%                               constructed by processSegmentEvents
%
%   Robel Yirdaw, December 2013
%       Modified, July 2014


    %If only one parent segment present, then it's a split event
    if (isnan(parentID(2)))
        %split
        %Set parent size
        parentSize = segSize;
        %Extract intensity amplitudes for children segments
        childInts = [segmentSetIN(childID(1),7);segmentSetIN(childID(2),7)];
        %Identify relative sizes for both
        biggerChildID = childID((childInts == max(childInts)));
        smallerChildID = childID(childID ~= biggerChildID);
        %We will use the parent's and the bigger child's size to determine
        %size probabilistically since the smaller one must be of size 1
        segIntensity = [segmentSetIN(parentID(1),7);segmentSetIN(biggerChildID(1),7)];
        
        %Call function to determine size for parent
        newParentSize = getSizeFromProb(parentSize,segIntensity,intensityInfo);
        
        %Set probabilistically determined size and calculate bigger child's
        %from the parent's size
        segmentSetOUT(parentID(1),8) = newParentSize;
        segmentSetOUT(biggerChildID(1),8) = newParentSize - 1;
        segmentSetOUT(smallerChildID(1),8) = 1;
    else
        %This is a merge event
        %Set child's size
        childSize_rev = segSize;        
        %Extract intensity amplitudes for parent segments        
        parentInts = [segmentSetIN(parentID(1),7);segmentSetIN(parentID(2),7)];                            
        %Identify relative sizes for both
        biggerParentID = parentID((parentInts == max(parentInts)));
        smallerParentID = parentID(parentID ~= biggerParentID(1));
        %We will use the child's and the bigger parent's size to determine
        %size probabilistically since the smaller one must be of size 1        
        segIntensity = [segmentSetIN(childID(1),7);segmentSetIN(biggerParentID(1),7)];
        
        %Call function to determine size for child
        newChildSize = getSizeFromProb(childSize_rev,segIntensity,intensityInfo);

        %Set probabilistically determined size and calculate bigger
        %parent's size from the child's size        
        segmentSetOUT(childID(1),8) = newChildSize;
        segmentSetOUT(biggerParentID(1),8) = newChildSize - 1;
        segmentSetOUT(smallerParentID(1),8) = 1;    
    end
    

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newSegSize = getSizeFromProb(segSize,segIntensity,intensityInfo)
%GETSIZEFROMPROB uses segments' intensity amplitude, the currently
%determined sizes and the intensity quantum to search for the most likely 
%size, via calcSegSize. This approach is used for cases where
%determining sizes based on events is not possible.
%   
%   INPUT:
%           segSize:            size of parent or child determined but to
%                               be adjusted here possible
%           segIntensity:       intensity amplitudes for 2 of the 3
%                               segments, depending on the event, as
%                               described in assignSizeFromProb
%           intensityInfo:      intensity quantum as defined in main sim.
%
%   OUTPUT:
%           newSegSize:         new size for segSize (possible that it may
%                               stay the same)
%
%   Robel Yirdaw, December 2013
%       Modified, July 2014


    %Initialize vars
    newSegSize = segSize;
    firstLowerProb = 0;
    searchDone = 0;
    
    %Calculate the probability of the segment having the size that was
    %determined based on its intensity and that of the second segment
    currProb = calcSegSizeProb(segSize(1),segIntensity(1),intensityInfo) *...
            calcSegSizeProb(segSize(1) - 1,segIntensity(2),intensityInfo);    
        
    %Next calcuate the immediate lower and higher probabilities possible,
    %by adjusting the size
    if (segSize > 2)
        firstLowerProb = calcSegSizeProb(segSize(1) - 1,segIntensity(1),intensityInfo) *...
            calcSegSizeProb(segSize(1) - 2,segIntensity(2),intensityInfo);
    end
    firstHigherProb = calcSegSizeProb(segSize(1) + 1,segIntensity(1),intensityInfo) *...
        calcSegSizeProb(segSize(1),segIntensity(2),intensityInfo);
    
    %Will perform search if the current probability either of the lower or
    %higher probabilities are greater. Otherwise, we are done
    if (firstLowerProb > currProb)
        currProb = firstLowerProb;
        sizeInc = -1;
        newSegSize = segSize(1) - 1;
    elseif (firstHigherProb > currProb)
        currProb = firstHigherProb;
        sizeInc = 1;
        newSegSize = segSize(1) + 1;
    else
        searchDone = 1;
    end

    %Perform search by repeatedly adjusting the size and recalculating the
    %probability of having that size
    while (~searchDone && (newSegSize > 1))
        %
        newProb = calcSegSizeProb(newSegSize + sizeInc,segIntensity(1),intensityInfo) *...
        calcSegSizeProb(newSegSize + sizeInc - 1,segIntensity(2),intensityInfo);
    
        if (newProb < currProb)
            searchDone = 1;
        else
            newSegSize = newSegSize + sizeInc;
        end
    end

end



function segSizeProb = calcSegSizeProb(segSize,segIntensity,intensityQuant)
%CALCSEGSIZEFROMPROB uses segments' intensity amplitude, the currently
%determined sizes and the intensity quantum to calculate the probability of
%having that size. This approach is used for cases where determining sizes 
%based on events is not possible.
%   
%   INPUT:
%           segSize:            size for which probabilit is to be
%                               calculated
%           segIntensity:       intensity amplitudes for 2 of the 3
%                               segments, depending on the event, as
%                               described in assignSizeFromProb
%           intensityInfo:      intensity quantum as defined in main sim.
%
%   OUTPUT:
%           segSizeProb:        probability for having the given size
%
%   Robel Yirdaw, December 2013
%       Modified, July 2014


    %meanIntensitySeg = segSize*intensityQuant(1);
    if (segIntensity > (segSize*intensityQuant(1)))
        segSizeProb = 1 - normcdf(segIntensity,segSize*intensityQuant(1),...
            sqrt(segSize)*intensityQuant(2));
    else        
        segSizeProb = normcdf(segIntensity,segSize*intensityQuant(1),...
            sqrt(segSize)*intensityQuant(2));
    end        

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function notRepIntSizeErr = perCaseRepairResult(tempRepFlag,parentID,...
    childID,segmentSetOUT,notRepIntSizeErr,caseStr,verbose)
%PERCASEREPAIRRESULT prints the repair outcome onto the screen and also
%returns the saves the string in notRepIntSizeErr for return.
%   
%   INPUT:
%           tempRepFlag:        a boolean indicating repair outcome
%           parentID:           segment ID of parents (1 or 2 element)
%           segmetSetOUT:       current table of segments and events
%                               constructed by processSegmentEvents
%           notRepIntSizeErr:   5 element cell with information for
%                               segements whose size was not repaired, as
%                               detailed in aggregStateFromCompTracks_new
%           caseStr:            a string specific corresponding to the 
%                               event case type
%           verbose:            boolean value indicating whether repair
%                               outcome will be printed to screen
%
%   OUTPUT:
%           notRepIntSizeErr:   updated 5 element cell with information for
%                               segements whose size was not repaired, as
%                               detailed in aggregStateFromCompTracks_new
%
%   Robel Yirdaw, December 2013
%       Modified, July 2014

    if (isnan(parentID(2)))
        %split
        if (tempRepFlag)
            outStr = sprintf('Repaired (%s): (%d:%d) -> (%d:%d) + (%d:%d). ',...
                caseStr,parentID(1),segmentSetOUT(parentID(1),8),...
                childID(1),segmentSetOUT(childID(1),8),childID(2),segmentSetOUT(childID(2),8));                                     
        elseif (~tempRepFlag)
            outStr = sprintf('Not repaired (%s): (%d:%d) -> (%d:%d) + (%d:%d). ',...
                caseStr,parentID(1),segmentSetOUT(parentID(1),8),childID(1),...
                segmentSetOUT(childID(1),8),childID(2),segmentSetOUT(childID(2),8));
            
            errMagnitude = segmentSetOUT(parentID(1),8) - ...
                    (segmentSetOUT(childID(1),8) + segmentSetOUT(childID(2),8));                           
                
            notRepIntSizeErr = setNotRepErrList(notRepIntSizeErr,parentID(1),...
                caseStr,outStr,'split',errMagnitude);
            
        end                                   
    else
        %merge
        if (tempRepFlag)                        
            outStr = sprintf('Repaired (%s): (%d:%d) + (%d:%d) -> (%d:%d). ',caseStr,...
                parentID(1),segmentSetOUT(parentID(1),8),parentID(2),...
                segmentSetOUT(parentID(2),8),childID(1),segmentSetOUT(childID(1),8));
        else
            outStr = sprintf('Not repaired (%s): (%d:%d) + (%d:%d) -> (%d:%d). ',...
                caseStr,parentID(1),segmentSetOUT(parentID(1),8),parentID(2),...
                segmentSetOUT(parentID(2),8),childID(1),segmentSetOUT(childID(1),8));
            
            errMagnitude = (segmentSetOUT(parentID(1),8) + ...
                    segmentSetOUT(parentID(2),8)) - segmentSetOUT(childID(1),8);
                
            notRepIntSizeErr = setNotRepErrList(notRepIntSizeErr,parentID,caseStr,...
                outStr,'merge',errMagnitude);
            
            fprintf(outStr);
        end
    end
    
    if (verbose)
        fprintf(outStr);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function notRepIntSizeErr = setNotRepErrList(notRepIntSizeErr,segID,errTypeStr,...
    outStr,eventType,errMagnitude)
%SETNOTREPERRLIST is used to populate notRepIntSizeErr.
%   
%   INPUT:
%           notRepIntSizeErr:   5 element cell with information for
%                               segements whose size was not repaired, as
%                               detailed in aggregStateFromCompTracks_new
%           segID:              ID of current segment
%           errTypeStr:         a string specific corresponding to the 
%                               event case type
%           outStr:             the string for repair outcome constructed
%                               in perCaseRepairResult
%           eventType:          a string indicating type of event, split or
%                               merge
%           errorMagnitude:     the amount of size discrepancy determined
%
%   OUTPUT:
%           notRepIntSizeErr:   updated 5 element cell with information for
%                               segements whose size was not repaired, as
%                               detailed in aggregStateFromCompTracks_new
%
%   Robel Yirdaw, December 2013
%       Modified, July 2014

    %First check if the segment is already in the list
    if (length(segID) == 2)
        segExists = (any(notRepIntSizeErr{1} == segID(1)) | ...
            any(notRepIntSizeErr{1} == segID(2)));
    else
        segExists = (any(notRepIntSizeErr{1} == segID(1)));
    end
    
    %The segment is not already in the list
    if (~segExists)
        tempIndx = sum(~isnan(notRepIntSizeErr{1})) + 1;
        notRepIntSizeErr{1}(tempIndx) = segID(1);
        notRepIntSizeErr{2}{tempIndx} = errTypeStr;
        notRepIntSizeErr{3}{tempIndx} = outStr;
        notRepIntSizeErr{4}{tempIndx,1} = eventType;
        notRepIntSizeErr{5}(tempIndx) = errMagnitude;
    end


end

