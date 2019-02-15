function plotTransportMotionAnalysisRes(motionAnalysisRes,indx2plot)

%extract toward and away fields from structure
structAway = motionAnalysisRes.awayFromRefPoint;
structToward = motionAnalysisRes.towardRefPoint; %#ok<NASGU>

%check whether there is a pause field
pauseExists = isfield(motionAnalysisRes,'pause');

%get names of supplied motion characteristics
motionChar = fieldnames(structAway);
numMotionChar = length(motionChar);

%% Plotting

for iMC = 1 : numMotionChar
    
    %get the instant speed distributions and statistics
    eval(['distrAway = structAway.' motionChar{iMC} '(indx2plot).values;'])
    eval(['statsAway = structAway.' motionChar{iMC} '(indx2plot).stats;'])
    eval(['distrToward = structToward.' motionChar{iMC} '(indx2plot).values;'])
    eval(['statsToward = structToward.' motionChar{iMC} '(indx2plot).stats;'])
    distrAll = [distrAway; distrToward];
    if strcmp(motionChar{iMC},'runTime') && pauseExists
        distrPause = motionAnalysisRes.pause.pauseTime(indx2plot).values;
        statsPause = motionAnalysisRes.pause.pauseTime(indx2plot).stats;
        distrAll = [distrAll; distrPause]; %#ok<AGROW>
    end
    
    %use the function "optimalHistogram" to get an estimate of the number of bins to
    %use in making the histogram
    n1 = optimalHistogram(distrAway,[],0);
    n2 = optimalHistogram(distrToward,[],0);
    numBins = mean([length(n1) length(n2)]);
    
    %find the minimum and maximum values and determine the bin locations
    minVal = min(distrAll);
    maxVal = max(distrAll);
    binDelta = (maxVal-minVal)/numBins;
    binPos = minVal:binDelta:maxVal;
    
    %now plot the histogram, mean, std and median
    n1 = hist(distrAway,binPos);
    n2 = hist(distrToward,binPos);
    n = max([n1 n2]);
    figure('Name',[motionChar{iMC} '  ' num2str(indx2plot)])
    if strcmp(motionChar{iMC},'runTime') && pauseExists
        numSubPlots = 3;
    else
        numSubPlots = 2;
    end
    subplot(1,numSubPlots,1)
    hist(distrAway,binPos)
    hold on
    plot(statsAway(1)*[1 1],[0 n+1],'g','LineWidth',2) %mean
    plot((statsAway(1)+statsAway(3))*[1 1],[0 n+1],'c--','LineWidth',2) %mean+1std
    plot((statsAway(1)-statsAway(3))*[1 1],[0 n+1],'c--','LineWidth',2) %mean-1std
    plot(statsAway(2)*[1 1],[0 n+1],'r','LineWidth',2) %median
    legend('distribution','mean','mean+1std','mean-1std','median')
    title('Away from reference point')
    subplot(1,numSubPlots,2)
    hist(distrToward,binPos)
    hold on
    plot(statsToward(1)*[1 1],[0 n+1],'g','LineWidth',2) %mean
    plot((statsToward(1)+statsToward(3))*[1 1],[0 n+1],'c--','LineWidth',2) %mean+1std
    plot((statsToward(1)-statsToward(3))*[1 1],[0 n+1],'c--','LineWidth',2) %mean-1std
    plot(statsToward(2)*[1 1],[0 n+1],'r','LineWidth',2) %median
    legend('distribution','mean','mean+1std','mean-1std','median')
    title('Toward reference point')
    if strcmp(motionChar{iMC},'runTime') && pauseExists
        subplot(1,numSubPlots,3)
        hist(distrPause,binPos)
        hold on
        plot(statsPause(1)*[1 1],[0 n+1],'g','LineWidth',2) %mean
        plot((statsPause(1)+statsPause(3))*[1 1],[0 n+1],'c--','LineWidth',2) %mean+1std
        plot((statsPause(1)-statsPause(3))*[1 1],[0 n+1],'c--','LineWidth',2) %mean-1std
        plot(statsPause(2)*[1 1],[0 n+1],'r','LineWidth',2) %median
        legend('distribution','mean','mean+1std','mean-1std','median')
        title('Pause time')
    end
    
end
