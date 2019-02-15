function seqOfEvents = removeSplitMergeArtifacts(seqOfEventsIn,replaceSegNum)

%copy sequence of events
seqOfEvents = seqOfEventsIn;

%find indices where there are splits
splitIndx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);

%go over all splits
while ~isempty(splitIndx)

    %get the first split
    iSplit = splitIndx(1);

    %find the two splitting segments and the time of splitting
    segment1 = seqOfEvents(iSplit,3);
    segment2 = seqOfEvents(iSplit,4);
    timeSplit = seqOfEvents(iSplit,1);

    %check whether these two segments merge with each other again
    iMerge = find(any(seqOfEvents(:,3:4)==segment1,2) & ...
        any(seqOfEvents(:,3:4)==segment2,2) & seqOfEvents(:,2)==2);
    
    %if they merge ...
    if ~isempty(iMerge)

        %calculate split to merge time
        timeMerge = seqOfEvents(iMerge,1);
        timeSplit2Merge = timeMerge - timeSplit;
        
        %if the split to merge time is only 1 frame
        if timeSplit2Merge == 1
            
            %discard both the split and the merge from the sequence of events
            seqOfEvents = seqOfEvents([1:iSplit-1 iSplit+1:iMerge-1 iMerge+1:end],:);
            
            %replace any mention of segment1 with segment2
            if replaceSegNum
                seqOfEvents(seqOfEvents(:,3)==segment1,3) = segment2;
                seqOfEvents(seqOfEvents(:,4)==segment1,4) = segment2;
            end
            
            %update indices of splits
            splitIndx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
            
        else
            
            %remove this split from list of splits
            splitIndx = splitIndx(2:end);
            
        end %(if timeSplit2Merge == 1 ... else ...)
        
    else
        
        %find when each of the segments ends
        iEnd1 = find(seqOfEvents(:,3)==segment1 & ...
            isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);
        timeEnd1 = seqOfEvents(iEnd1,1);
        iEnd2 = find(seqOfEvents(:,3)==segment2 & ...
            isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);
        timeEnd2 = seqOfEvents(iEnd2,1);
        
        %get the split-to-end time for both segments
        timeSplit2End1 = timeEnd1 - timeSplit;
        timeSplit2End2 = timeEnd2 - timeSplit;
        
        %if either time is equal to 0 frame, discard the split
        if any([timeSplit2End1 timeSplit2End2]==0)

            %if the splitting segment (segment1) ends immediately
            if timeSplit2End1 == 0

                %discard both the split and the end of segment1
                seqOfEvents = seqOfEvents([1:iSplit-1 iSplit+1:iEnd1-1 iEnd1+1:end],:);

            elseif timeSplit2End2 == 0 %if the segment that got split from (segment2) ends immediately
                
                %discard both the split and the end of segment2
                seqOfEvents = seqOfEvents([1:iSplit-1 iSplit+1:iEnd2-1 iEnd2+1:end],:);
                
                %replace any mention of segment1 with segment2
                if replaceSegNum
                    seqOfEvents(seqOfEvents(:,3)==segment1,3) = segment2;
                    seqOfEvents(seqOfEvents(:,4)==segment1,4) = segment2;
                end

            end

            %update indices of splits
            splitIndx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);

        else

            %remove this split from list of splits
            splitIndx = splitIndx(2:end);

        end

    end %(if ~isempty(iMerge) ... else ...)
    
end %(while ~isempty(splitIndx))

%find indices where there are merges
mergeIndx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);

%go over all merges
while ~isempty(mergeIndx)

    %get the first merge
    iMerge = mergeIndx(1);

    %find the two merging segments and the time of merging
    segment1 = seqOfEvents(iMerge,3);
    segment2 = seqOfEvents(iMerge,4);
    timeMerge = seqOfEvents(iMerge,1);

    %find when each of the segments started
    iStart1 = find(seqOfEvents(:,3)==segment1 & ...
        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
    timeStart1 = seqOfEvents(iStart1,1);
    iStart2 = find(seqOfEvents(:,3)==segment2 & ...
        isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==1);
    timeStart2 = seqOfEvents(iStart2,1);

    %get the start-to-merge time for both segments
    timeStart12Merge = timeMerge - timeStart1;
    timeStart22Merge = timeMerge - timeStart2;

    %if either time is equal to 1 frame, discard the merge
    if any([timeStart12Merge timeStart22Merge]==1)

        %if the merging segment (segment1) started just before
        if timeStart12Merge == 1

            %discard both the merge and the start of segment1
            seqOfEvents = seqOfEvents([1:iStart1-1 iStart1+1:iMerge-1 iMerge+1:end],:);

        elseif timeStart22Merge == 1 %if the segment that got merged with (segment2) started just before

            %discard both the merge and the start of segment2
            seqOfEvents = seqOfEvents([1:iStart2-1 iStart2+1:iMerge-1 iMerge+1:end],:);

            %replace any mention of segment1 with segment2
            if replaceSegNum
                seqOfEvents(seqOfEvents(:,3)==segment1,3) = segment2;
                seqOfEvents(seqOfEvents(:,4)==segment1,4) = segment2;
            end

        end

        %update indices of merges
        mergeIndx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2)==2);

    else

        %remove this merge from list of merges
        mergeIndx = mergeIndx(2:end);

    end

end %(while ~isempty(mergeIndx))

