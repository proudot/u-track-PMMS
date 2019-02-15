function mismatchInfo_sims = plotAggregStateMismatchInfo(mismatchInfo)
%PLOTAGGREGSTATEMISMATCHINFO takes a set of segment size mismatch
%information from multiple simulations and produces a figure summarizing
%the mismatches, over the set of simulations. %   A single figure with 5 
%subplots will be created.
%   
%
%   INPUT:
%       mismatchInfo:       a struct array with mismatchInfo for each
%                           simulation. mismatchInfo is a struct
%                           created by checkAggregState.
%
%       mismatchInfo_sims:  a struct of summary quantities calculated
%                           and plotted here. Fields are:
%                           1) numSims: number of simulations processed
%                           2) numSimsWithMismatches: number of simulations
%                           that contained mismatches
%                           3) numIters: number of iterations that 
%                           contained mismatches
%                           4) largestClustSize: largest cluster for all
%                           simulations
%                           5) minDev:  the smallest size deviation
%                           6) maxDev:  the largest size deviation
%                           7) totMismatches: total number of mismatches
%                           over all simulations
%                           8) diffNumClust:  difference in number of
%                           clusters at each iteration
%                           9) clustDist_AS:  number of clusters, by
%                           cluster size, at each iteration from
%                           agggregState
%                           10) clustDist_AS_count:  number of iterations at
%                           which the sum in #8 was performed
%                           11) clustDist_R2C:  number of clusters, by
%                           cluster size, at each iteration from
%                           recept2clustAssign
%                           12) clustDist_R2C_count:  number of iterations at
%                           which the sum in #10 was performed
%                           13) meanClustDist_AS:  mean number of cluster
%                           from aggregState
%                           14) meanClustDist_R2C:  mean number of cluster
%                           from aggregState
%                           15) diffPerClustSize_Hist: histogram of
%                            mismatches by cluster size
%                           16) diffPerClustSize_Hist: histogram of size
%                           differences as a fraction of total mismatches
%                           17) bSeriesVals: values plotted on bar series
%                           18) bSeriesXVals: x values plotted on bar
%                           series which are mismatch amounts
%                           19) totDevsBySize: total number of mismatches
%                           by cluster size
%                           20) totDevsBySize_frac: total number of 
%                           mismatches by cluster size as a fraction
%                           21) mismatchDiffsPerIter: mismatch amounts per
%                           iterations
%                           22) mismatchDiffsPerIter_frac: mismatch amounts per
%                           per iterations as fractions
%                           23) diffPerClustSize_Quant: difference in sizes
%                           by cluster size at each iteration
%                           24) diffPerClustSize_Quant_norm: difference in 
%                           sizes by cluster size at each iteration
%                           normalized by true number of cluster from R2C
%                           25) diffPerClustSize_Quant_absSum: sum of the 
%                           absoulte value of differences by cluster size
%
%   Robel Yirdaw, January 2014
%



    mismatchInfo_sims = [];
    
    if (all(isnan([mismatchInfo.totMismatches])))
        fprintf('\n==================plotAggregStateMismatchInfo==================');
        fprintf('\nNo mismatches found in mismatchInfo.');
        fprintf('\n===============================================================');
    else
        %Determine number of simulations
        numSims = length(mismatchInfo);        
        %The following used for plots 2b and 5.
        numSimWithMismatches = sum(~isnan([mismatchInfo.totMismatches]));
        numIters = length([mismatchInfo.diffNumClust]);
        simLargestClustSize = max([mismatchInfo.largestClustSize]);       
        %Minimum, maximum and largest absolute deviations
        minDev = min([mismatchInfo.minDev]);
        maxDev = max([mismatchInfo.maxDev]);
        largestAbsDev = max(abs(minDev),maxDev);
        
        %Total number of mismatches
        %Must use nansum since some in the array of mismatchInfo entries
        %may not have mismatches (multi sim case).
        totMismatches = nansum([mismatchInfo.totMismatches]);
        
        %Intialize variables
        diffNumClust = NaN(numIters,1);
        clustDist_AS = NaN(simLargestClustSize,numIters);
        clustDist_AS_count = NaN(simLargestClustSize,1);
        clustDist_R2C = NaN(simLargestClustSize,numIters);
        clustDist_R2C_count = NaN(simLargestClustSize,1);
        diffPerClustSize_Hist = NaN(max(abs(minDev),maxDev),2);
        %Will be plotting a bar series for subplot 1
        bSeriesVals = NaN(numIters,abs(minDev)+maxDev+1);
        
        totDevsBySize = NaN(simLargestClustSize,1);
        totDiffPerClustSize_Quant = NaN(simLargestClustSize,2);    
        diffPerClustSize_Quant_norm = NaN(simLargestClustSize,numIters,numSims);
        diffPerClustSize_Quant_absSum = NaN(simLargestClustSize,numIters);

        mismatchDiffsPerIter = NaN(largestAbsDev,1);

        for simCount=1:numSims
            %Since some of the sims. can be mismatch free, do the following
            %only for those with mismatches.
            if (~isnan(mismatchInfo(simCount).totMismatches))
                %Get the current largest cluster size
                currLargestClustSize = mismatchInfo(simCount).largestClustSize;
                %Difference in number of clusters at each iteration will be
                %plotted using a barseries
                diffNumClust = nansum([diffNumClust mismatchInfo(simCount).diffNumClust],2);
                if (simCount == 1)
                    bSeriesVals = histc(mismatchInfo(simCount).diffPerClustSize_Quant,minDev:maxDev)';
                else
                    bSeriesVals = bSeriesVals +...
                        histc(mismatchInfo(simCount).diffPerClustSize_Quant,minDev:maxDev)';
                end

                currHistLength = [length(mismatchInfo(simCount).diffPerClustSize_Hist(:,1))...
                    length(mismatchInfo(simCount).diffPerClustSize_Hist(:,2))];
                
                %Accumulate difference per cluster size
                diffPerClustSize_Hist(1:currHistLength(1),1:2) =...
                    [nansum([diffPerClustSize_Hist(1:currHistLength(1),1)...
                        mismatchInfo(simCount).diffPerClustSize_Hist(1:currHistLength(1),1)],2) ...
                    nansum([diffPerClustSize_Hist(1:currHistLength(2),2)...
                        mismatchInfo(simCount).diffPerClustSize_Hist(1:currHistLength(2),2)],2)];
                    
                %Accumulate total devaitions by cluster size    
                totDevsBySize(1:currLargestClustSize,1) =...
                    nansum([totDevsBySize(1:currLargestClustSize,1)...
                    nansum(mismatchInfo(simCount).diffPerClustSize_Quant ~= 0 &...
                        ~isnan(mismatchInfo(simCount).diffPerClustSize_Quant),2)],2);
                    
                %Mismatch differences per iteration
                currDiffsLength = length(mismatchInfo(simCount).mismatchDiffsPerIter);
                mismatchDiffsPerIter(1:currDiffsLength,1) =...
                    nansum([mismatchDiffsPerIter(1:currDiffsLength,1)...
                    mismatchInfo(simCount).mismatchDiffsPerIter],2);
                                
                for clustSize=1:simLargestClustSize
                    if (clustSize <= mismatchInfo(simCount).largestClustSize)
                        %Accumulate absolute value of differences
                        diffPerClustSize_Quant_absSum(clustSize,:) =... 
                            nansum([diffPerClustSize_Quant_absSum(clustSize,:);... 
                            abs(mismatchInfo(simCount).diffPerClustSize_Quant(clustSize,:))]);
                        
                        %Calculate normalized differences per cluster size
                        diffPerClustSize_Quant_norm(clustSize,1:numIters,simCount) =...                        
                            mismatchInfo(simCount).diffPerClustSize_Quant(clustSize,:)./...
                            mismatchInfo(simCount).clustDist_R2C(clustSize,:);
                        %Number of clustrs, by cluster size, at each
                        %iteration from aggregState and R2C
                        clustDist_AS(clustSize,:) = nansum([clustDist_AS(clustSize,:);...
                            mismatchInfo(simCount).clustDist_AS(clustSize,:)]);
                        clustDist_R2C(clustSize,:) = nansum([clustDist_R2C(clustSize,:);...
                            mismatchInfo(simCount).clustDist_R2C(clustSize,:)]);  
                        
                        %Accumulate the number of iterations with clustSize
                        %i.e. number of data points for the above two sums
                        clustDist_AS_count(clustSize,1) = nansum([clustDist_AS_count(clustSize,1)...
                            sum(~isnan(mismatchInfo(simCount).clustDist_AS(clustSize,:)),2)],2);
                        clustDist_R2C_count(clustSize,1) = nansum([clustDist_R2C_count(clustSize,1)...
                            sum(~isnan(mismatchInfo(simCount).clustDist_R2C(clustSize,:)),2)],2);                        
                    end

                end

            end %if current sim has mismatches
            
        end %for simCount

        for sizeIndx=1:simLargestClustSize
            totDiffPerClustSize_Quant(sizeIndx,1) =...
                nanmean(abs(diffPerClustSize_Quant_norm(sizeIndx,...
                ~isinf(diffPerClustSize_Quant_norm(sizeIndx,:,:)) ) ) );
        end            
        totDiffPerClustSize_Quant(:,2) = nansum(diffPerClustSize_Quant_absSum,2)./...
            nansum(clustDist_R2C,2);     


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        figure();
        subplot(3,2,1:2);
        [AX,H1,H2] = plotyy(1:numIters,bSeriesVals,1:numIters,diffNumClust,'bar','plot');
        set(H1,'BarLayout','stacked');
        set(get(AX(1),'Ylabel'),'String','Mismatch Magnitudes');
        set(get(AX(2),'Ylabel'),'String','Cluster # Mismatches');    
        set(H2,'Marker','o','Color','b');      
        xlabel(gca,'Iteration');
        legendVals = minDev:maxDev;
        for dataCount=1:length(legendVals)
            set(H1(dataCount),'DisplayName',num2str(legendVals(dataCount)));
        end
        set(H2(1),'DisplayName','Cluster #');
        lH = legend(AX(1),'Location','East','Orientation','Verictal');
        set(lH,'FontSize',8);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(3,2,3);

        meanClustDist_AS = nansum(clustDist_AS,2)./clustDist_AS_count;
        meanClustDist_R2C = nansum(clustDist_R2C,2)./clustDist_R2C_count;

        [AX,H1,H2] = plotyy(1:length(totDiffPerClustSize_Quant(:,1)),totDiffPerClustSize_Quant,...
            1:length(meanClustDist_AS),[meanClustDist_R2C meanClustDist_AS],'bar','plot');
        set(get(AX(1),'Ylabel'),'String','(rel. error)');
        set(get(AX(2),'Ylabel'),'String','mean (count)');    
        xlabel(gca,'Cluster Size');  
        set(H1(1),'FaceColor',[0;0;0])
        set(H2(1),'Marker','sq','Color','k','LineStyle',':');
        set(H2(2),'Marker','+','Color','b','LineStyle',':');
        set(H1(1),'FaceColor',[1;0.6275;0.4784])
        set(H1(2),'FaceColor',[0.4784;1;0.9216])

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(3,2,4);
        
        totDevsBySize_frac = totDevsBySize./totMismatches;
        bar(1:length(totDevsBySize_frac),totDevsBySize_frac.*100);
        xlabel(gca,'Cluster Size');
        ylabel(gca,'% (total mismatches)');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        subplot(3,2,5);
        
        diffPerClustSize_Hist_frac(:,1) = [(-1:-1:minDev)';(1:maxDev)'];
        diffPerClustSize_Hist_frac(:,2) = ([diffPerClustSize_Hist(1:abs(minDev),1);...
            diffPerClustSize_Hist(1:maxDev,2)])./totMismatches;
        bar(diffPerClustSize_Hist_frac(:,1),...
            diffPerClustSize_Hist_frac(:,2).*100);    
        xlabel(gca,'Deviation');  
        ylabel(gca,'% (total mismatches)');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        
        subplot(3,2,6);
        mismatchDiffsPerIter_frac = mismatchDiffsPerIter./(numIters*numSimWithMismatches);
        bar(1:max(abs(minDev),maxDev),mismatchDiffsPerIter_frac.*100);
        xlabel(gca,'abs(Deviation)');    
        ylabel(gca,'% (iterations)');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
        %Save calculated values for return
        mismatchInfo_sims.numSims = numSims;
        mismatchInfo_sims.numSimWithMismatches = numSimWithMismatches;
        mismatchInfo_sims.numIters = numIters;	
        mismatchInfo_sims.largestClustSize = simLargestClustSize;
        mismatchInfo_sims.minDev = minDev;
        mismatchInfo_sims.maxDev = maxDev;
        mismatchInfo_sims.totMismatches = totMismatches;
        mismatchInfo_sims.diffNumClust = diffNumClust;
        mismatchInfo_sims.clustDist_AS = clustDist_AS;
        mismatchInfo_sims.clustDist_AS_count = clustDist_AS_count;
        mismatchInfo_sims.clustDist_R2C = clustDist_R2C;
        mismatchInfo_sims.clustDist_R2C_count = clustDist_R2C_count;
        mismatchInfo_sims.meanClustDist_AS = meanClustDist_AS;
        mismatchInfo_sims.meanClustDist_R2C = meanClustDist_R2C;
        mismatchInfo_sims.diffPerClustSize_Hist = diffPerClustSize_Hist;
        mismatchInfo_sims.diffPerClustSize_Hist_frac = diffPerClustSize_Hist_frac;
        mismatchInfo_sims.bSeriesVals = bSeriesVals;
        mismatchInfo_sims.bSeriesXvals = legendVals;
        mismatchInfo_sims.totDevsBySize = totDevsBySize;
        mismatchInfo_sims.totDevsBySize_frac = totDevsBySize_frac;
        mismatchInfo_sims.mismatchDiffsPerIter = mismatchDiffsPerIter;
        mismatchInfo_sims.mismatchDiffsPerIter_frac = mismatchDiffsPerIter_frac;
        mismatchInfo_sims.totDiffPerClustSize_Quant = totDiffPerClustSize_Quant;
        mismatchInfo_sims.diffPerClustSize_Quant_norm = diffPerClustSize_Quant_norm;
        mismatchInfo_sims.diffPerClustSize_Quant_absSum = diffPerClustSize_Quant_absSum;

    end %if mismatches exist
    
end
