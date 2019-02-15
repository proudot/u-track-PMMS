function [percentageGoodLink,percentageBadLink,percentageMissingLink,noise_var_MSE, schemesStat,cpuTime,tracksFinal] ... 
        = track_stat(movieInfo,tracksSim,trackDiffCoef,costMatrices,kalmanFunctions,gapCloseParam,saveFolder,movieParam,schemeName)
    
    probDim=2;
    verbose=1;
    saveResults.dir =  [saveFolder schemeName]; %directory where to save input and output
    mkdir(saveResults.dir);
    
    watch_KF_iter=0;
    
% $$$     revMovieParam=movieParam
% $$$     revMovieParam.imageDir=[movieParam.imageDir '/rev/']
% $$$     
    starttime=cputime;    
    [tracksFinal,kalmanInfoLink,errFlag] = ...
        trackCloseGapsKalmanSparse(movieInfo,costMatrices, ... 
                                     gapCloseParam,kalmanFunctions,...
                                     probDim,saveResults,verbose,'iter_link_nb',1,   ...
                                     'watch_KF_iter',watch_KF_iter,'movieParam',movieParam);
    
    

    % computational time
    cpuTime=cputime-starttime;

    %linking percentage
    linkStats=scoreLinksGapsMSLinkType(tracksFinal,tracksSim);
    avgperc= nanmean(linkStats(:,3,:)./linkStats(:,1,:))*100;
    percentageGoodLink=avgperc;
    avgperc= nanmean(linkStats(:,4,:)./linkStats(:,1,:))*100;
    percentageBadLink=avgperc;
    avgperc= nanmean( (linkStats(:,1,:)-linkStats(:,3,:))./linkStats(:,1,:))*100;
    percentageMissingLink=avgperc;
    
    %noise variance estimate
    % !!! does not work with merge and split.
    a_noise_var_MSE=0;
    timeStep=0.1;
    stepStd = sqrt(2*trackDiffCoef*timeStep);
    for track_idx= 1:size(tracksSim,1)
        process_noise_GT=stepStd(track_idx);
        if ndims(kalmanInfoLink(1).noiseVar)>3
            noise_vars=getTrackIndeNoiseVars(kalmanInfoLink,tracksSim,track_idx,3);
        else
            noise_vars=getTrackNoiseVars(kalmanInfoLink,tracksSim,track_idx);
        end 
        a_noise_var_MSE=a_noise_var_MSE+sum((squeeze(noise_vars(1,1,:)-process_noise_GT).^2));
    end 
    noise_var_MSE=a_noise_var_MSE/size(tracksSim,1);
    
    
    % Schemes stats
    schemes=[];
    for track_idx= 1:size(tracksFinal,1)
        schemes=[schemes; getTrackBestSchemes(kalmanInfoLink,tracksFinal,track_idx)];
    end
    minScheme=1;
    maxScheme=2;

    schemesStat=nan(maxScheme-minScheme+1,1);
    for s = minScheme:maxScheme
        schemesStat(s)=sum(schemes==s);
    end
end