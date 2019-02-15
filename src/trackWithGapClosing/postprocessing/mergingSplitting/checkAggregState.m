function [mismatchInfo,aggregStateFromR2C,aggregStateALL,compTracksFromR2C] =...
    checkAggregState(compTracksALT,r2c,simIndx)
%CHECKAGGREGSTATE compares sizes in aggregState which are determined from
%tracks to sizes that are determined directly from recept2clustAssign, the
%ground truth.
%
%   INPUT:
%           compTracksALT:      compTracks in the alternative format
%           r2c:                recept2clustAssign output by the simulation
%                               function, receptorAggregationSimple_new
%           simIndx:            simulation number
%
%   OUTPUT:
%           mismatchInfo:       a struct with size mismatch information
%                               with fields:
%                               1) misMatchIters: iterations with size
%                               mismatches detected
%                               2) diffNumClust:  difference in number of
%                               clusters present at each iteration
%                               3) diffPerClustID: difference in number of
%                               clusters present by cluster ID 
%                               4) diffPerClustSize_Quant: difference in
%                               number of clusters by cluster size at each
%                               iteration
%                               5) diffPerClustSize_ID: difference in
%                               number of clusters by cluster ID at each
%                               iteration (NOT USED)
%                               6) totDiffPerClustSize_Quant: mean of
%                               differences in number of cluster by cluster
%                               size
%                               7) diffPerClustSize_Hist: histogram of
%                               mismatches by cluster size
%                               8) clustDist_AS:  number of clusters, by
%                               cluster size, at each iteration from
%                               agggregState
%                               9) clustDist_R2C:  number of clusters, by
%                               cluster size, at each iteration from
%                               recept2clustAssign
%                               10) mismatchDiffsPerIter: mismatch amounts
%                               per iterations
%                               11) largestClustSize: the largest cluster
%                               size in this simulation
%                               12) minDev:  the smallest size deviation
%                               13) maxDev:  the largest size deviation
%                               14) totMismatches: total number of size
%                               mismatches detected
%
%   Robel Yirdaw, January 2014
%       Modified, 08/20/14
%



    fprintf('\n========================checkAggregState========================');
    fprintf('\nComparing sizes from aggregState and receptor2cluster.');
    %Determie number of receptors and simulation iterations
    [numReceptors,numIters] = size(r2c);
    
    %aggregState will be constructed from recept2clustAssign. aggregState
    %from compTracks will also be adjusted by removing blank rows and
    %returned (will match that from r2c).    
    aggregStateFromR2C = NaN(numReceptors,numIters);   
    aggregStateALL = NaN(numReceptors,numIters);
    
    %Initialize flags, counters and struct
    errorFlag = 0;
    mismatchCount = 0;
    mismatchInfo = struct('misMatchIters',NaN(numIters,2),'diffNumClust',NaN(numIters,1),...
        'diffPerClustID',NaN(numReceptors*2,numIters),'diffPerClustSize_Quant',NaN(numReceptors*2,numIters),...
        'diffPerClustSize_ID',NaN(numReceptors*2,numIters),...
        'totDiffPerClustSize_Quant', NaN(numReceptors*2,2),...
        'diffPerClustSize_Hist',NaN(numReceptors,2),...
        'clustDist_R2C',NaN(numReceptors*2,numIters),...
        'clustDist_AS',NaN(numReceptors*2,numIters),...
        'mismatchDiffsPerIter',NaN(numReceptors,1),...
        'largestClustSize',NaN,'minDev',NaN,'maxDev',NaN,'totMismatches',NaN);
    
    %The number of tracks present
    numTracks = length(compTracksALT.defaultFormatTracks);
    
    %Will construct a compTracks with aggregState from R2C - useful for
    %other functions that take compTracks
    compTracksFromR2C = compTracksALT;
    compTracksFromR2C.alternativeFormatTracks = [];
    for trackCount=1:numTracks
        [tempRows,tempCols] = size(compTracksFromR2C.defaultFormatTracks(trackCount).aggregState);
        compTracksFromR2C.defaultFormatTracks(trackCount).aggregState = NaN(tempRows,tempCols);        
        clear tempRows tempCols
    end
        
    %The check is performed at each iteration
    largestClustSize = 0;
    for iterCount=1:numIters
        %Largest cluster
        numClustersR2C = max(r2c(:,iterCount));
        %Tally size for each cluster (number of receptors) and construct
        %aggregStateFromR2C
        for clusterID=1:numClustersR2C
            aggregStateFromR2C(clusterID,iterCount) = sum(r2c(:,iterCount) == clusterID);
        end
        
        %Corresponding values in aggregState, from all tracks
        colFromAggregState = NaN(numReceptors,1);
        for tracksIter=1:numTracks
            tempAggregState = compTracksALT.defaultFormatTracks(tracksIter).aggregState;
            tempTracksFeatIndxCG = compTracksALT.defaultFormatTracks(tracksIter).tracksFeatIndxCG;
            
            %Get sizes at current iteration
            %Modified 082014 to allow sparse compTracks/aggregState wich
            %substitutes NaNs with 0s. This assumes that a segment will not
            %have size 0.
            currVals = tempAggregState(~isnan(tempAggregState(:,iterCount)) &...
                (tempAggregState(:,iterCount) ~= 0),iterCount);
            %Get indices from tracksFeatIndxCG
            currIndx = tempTracksFeatIndxCG(tempTracksFeatIndxCG(:,iterCount)~=0,iterCount);
            %Set up the column vector of sizes
            colFromAggregState(currIndx,1) = currVals;
            
            %Need to form a logical indices here
            aggregStateIndxLogical = (~isnan(tempAggregState(:,iterCount)) &...
                (tempAggregState(:,iterCount) ~= 0) );
            %Assign aggregStateFromR2C
            compTracksFromR2C.defaultFormatTracks(tracksIter).aggregState(aggregStateIndxLogical,iterCount) =...
                aggregStateFromR2C(currIndx,iterCount);
            
            clear currIndx currVals tempAggregState...
                tempTracksFeatIndxCG aggregStateIndx
        end
        %Number of clusters from aggregState
        numClustersAS = sum(~isnan(colFromAggregState));

        %Compare with values from r2c
        try       
            %diff = aggregStateFromR2C(1:tempMax,iterCount) - colFromAggregState;
            diff = sort(aggregStateFromR2C(1:numClustersR2C,iterCount)) - ...
                sort(colFromAggregState(1:numClustersAS,1));
            if (any(diff))
                %Save difference in number of clusters
                mismatchInfo.diffNumClust(iterCount,1) = numClustersAS - numClustersR2C;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Next determine difference in cluster sizes of individual
                %clusters/segments.
                
                tempMaxClustID = max(numClustersR2C,numClustersAS);
                for clustID=1:tempMaxClustID
                    tempR2Cval = 0;
                    if (~isnan(aggregStateFromR2C(clustID,iterCount)))
                        tempR2Cval = aggregStateFromR2C(clustID,iterCount);
                    end
                    tempASval = 0;
                    if (~isnan(colFromAggregState(clustID,1)))
                        tempASval = colFromAggregState(clustID,1);
                    end
                    tempDiff = tempASval - tempR2Cval; 
                    %Save difference by ID.
                    mismatchInfo.diffPerClustID(clustID,iterCount) = tempDiff;
                    
                    %Save by size also. The current cluster/seg has size
                    %tempR2Cval.
                    %mismatchInfo.diffPerClustSize_ID(tempR2Cval,iterCount) = tempDiff;            

                    clear tempDiff tempR2Cval tempASval tempCol
                end      
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Determine difference in cluster qunatities of all
                %cluster sizes.
                tempMaxClustSize = max(max(aggregStateFromR2C(1:numClustersR2C,iterCount)),...
                    max(colFromAggregState(1:numClustersAS,1)));     
                if (tempMaxClustSize > largestClustSize)
                    largestClustSize = tempMaxClustSize;
                end                
                for clustSizeIndx=1:tempMaxClustSize
                    if (any(aggregStateFromR2C(1:numClustersR2C,iterCount) == clustSizeIndx) ||...
                            any(colFromAggregState(1:numClustersAS,1) == clustSizeIndx) )
                        %The current cluster size is found in r2c
                        %Quantity:
                        numClustCurrSize_R2C =...
                            sum(aggregStateFromR2C(1:numClustersR2C,iterCount) == clustSizeIndx);
                        %Get quantity from aggregState to compare
                        numClustCurrSize_AS =...
                            sum(colFromAggregState(1:numClustersAS,1) == clustSizeIndx);
                        
                        diffCurrSize = numClustCurrSize_AS - numClustCurrSize_R2C;
                        mismatchInfo.diffPerClustSize_Quant(clustSizeIndx,iterCount) =...
                            diffCurrSize;
                        
                        mismatchInfo.clustDist_R2C(clustSizeIndx,iterCount) =...
                            numClustCurrSize_R2C;
                        mismatchInfo.clustDist_AS(clustSizeIndx,iterCount) =...
                            numClustCurrSize_AS;
                        
                        %Determine histogram of size differences
                        tempCol = 2;
                        if (diffCurrSize < 0)
                            tempCol = 1;
                        elseif (diffCurrSize > 0)
                            tempCol = 2;
                        end
                        if (diffCurrSize ~= 0)
                            if (isnan(mismatchInfo.diffPerClustSize_Hist(abs(diffCurrSize),tempCol)))
                                mismatchInfo.diffPerClustSize_Hist(abs(diffCurrSize),tempCol) = 1;
                            else
                                mismatchInfo.diffPerClustSize_Hist(abs(diffCurrSize),tempCol) = ...
                                    mismatchInfo.diffPerClustSize_Hist(abs(diffCurrSize),tempCol) + 1;
                            end
                        end                        

                        clear tempCol diffCurrSize numClustCurrSize_R2C numClustCurrSize_AS
                    end
                end %for
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (~errorFlag)
                    %Found mismatch in sizes - dump curr. values
                    errorFlag = 1;          
                    mismatchCount = mismatchCount + 1;
                    mismatchInfo.misMatchIters(mismatchCount,1) = iterCount;

                    %fprintf('\n==================aggregStateCompTracks test script==================');
                    fprintf('\nSize mismatch at iteration %d (tracksAmp col %d), sim. %d.',...
                        iterCount,(iterCount-1)*8 + 4,simIndx);
                    %{
                    fprintf('\nSizes from receptor2cluster: ');
                    aggregStateFromR2C_compact(1:numClustersR2C,iterCount)'
                    fprintf('Sizes from aggregState: ');
                    colFromAggregState'
                    %fprintf('\n===================================================================\n\n');                     
                    %pause;
                    %Uncomment line below to terminate on error
                    %return;                
                %}
                end
            elseif (errorFlag && ~any(diff))
                mismatchInfo.misMatchIters(mismatchCount,2) = iterCount;
                errorFlag = 0;
            end
        
        catch exception
            disp(exception);
            
            %The check above resulted in an error - dump curr. values
            errorFlag = 2;            
            
            %fprintf('\n==================aggregStateCompTracks test script==================');
            fprintf('\nSize check exception at iteration %d, sim. %d.\n',iterCount,simIndx);
            fprintf('Sizes from receptor2cluster: ');
            aggregStateFromR2C(1:numClustersR2C,iterCount)'
            fprintf('Sizes from aggregState: ');
            colFromAggregState'
            %fprintf('\n===================================================================\n\n');            
            pause;
            %Uncomment line below to terminate on error
            %return;

        end %try catch
        
        %save sorted aggregState
        aggregStateALL(1:numReceptors,iterCount) = colFromAggregState;
        
        
        clear colFromAggregState numClustersR2C numClustersAS tempMaxClustSize
    end %for each iteration

    mismatchInfo.misMatchIters(mismatchCount+1:end,:) = [];
    
    if (mismatchCount == 0)
        fprintf('\nNo size mismatches found (sim. %d).',simIndx);
    else      
        %Details on mismatches
        
        for sizeIndx=1:largestClustSize
            tempNormDiffs = mismatchInfo.diffPerClustSize_Quant(sizeIndx,:)./...
                mismatchInfo.clustDist_R2C(sizeIndx,:);
            mismatchInfo.totDiffPerClustSize_Quant(sizeIndx,1) =...
                nanmean(abs(tempNormDiffs(tempNormDiffs ~= inf)));    
            
            tempSumDiff = nansum(abs(mismatchInfo.diffPerClustSize_Quant(sizeIndx,:)));
            tempSumQuant = nansum(mismatchInfo.clustDist_R2C(sizeIndx,:));
            mismatchInfo.totDiffPerClustSize_Quant(sizeIndx,2) = tempSumDiff/tempSumQuant;
            
            %mismatchInfo.totDiffPerClustSize_Quant(sizeIndx,2) =...
            %    nanmean(tempNormVals);
            
            clear tempNormVals tempSumDiff tempSumQuant
        end
              
        
        mismatchInfo.largestClustSize = largestClustSize;
        mismatchInfo.minDev = min(min(mismatchInfo.diffPerClustSize_Quant));
        mismatchInfo.maxDev = max(max(mismatchInfo.diffPerClustSize_Quant));
        mismatchInfo.totMismatches =...
            sum(sum(~isnan(mismatchInfo.diffPerClustSize_Quant) &...
            (mismatchInfo.diffPerClustSize_Quant ~= 0))); 
        
        largestAbsDev = max(max(abs(mismatchInfo.minDev),mismatchInfo.maxDev));
        for diffIndx=1:largestAbsDev
            currDiffIters = any(abs(mismatchInfo.diffPerClustSize_Quant) == diffIndx);
            if (any(currDiffIters))
                mismatchInfo.mismatchDiffsPerIter(diffIndx,1) = sum(currDiffIters);
            end
            clear currDiffIters
        end

        %trim
        mismatchInfo.totDiffPerClustSize_Quant(largestClustSize+1:end,:) = [];
        %mismatchInfo.diffPerClustSize_ID(largestClustSize+1:end,:) = [];
        mismatchInfo.diffPerClustSize_Quant(largestClustSize+1:end,:) = [];
        mismatchInfo.clustDist_R2C(largestClustSize+1:end,:) = [];
        mismatchInfo.clustDist_AS(largestClustSize+1:end,:) = [];
        mismatchInfo.diffPerClustSize_Hist(...
            max(abs(mismatchInfo.minDev),mismatchInfo.maxDev)+1:end,:) = [];
        mismatchInfo.mismatchDiffsPerIter(largestAbsDev+1:end,:) = [];
        
    end
    
    fprintf('\n=====================================================================\n\n');                    

end %function

