function [trackStateNoiseVars] = getTrackNoiseVars (kalmanInfoLink,tracksFinal,trackIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;

trackStateNoiseVars=NaN(4,4,size(tracksFinal(trackIdx).tracksFeatIndxCG,2));
for t = startTime:endTime
    featureIdx=tracksFinal(trackIdx).tracksFeatIndxCG(tIdx);
    
    if featureIdx ~=0    
        trackStateNoiseVars(:,:,tIdx)=kalmanInfoLink(t).noiseVar(:, ...
                                                          :,tracksFinal(trackIdx).tracksFeatIndxCG(tIdx));

    else
        trackStateNoiseVars(:,:,tIdx)=zeros(4);
    end 
    tIdx=tIdx+1;
end