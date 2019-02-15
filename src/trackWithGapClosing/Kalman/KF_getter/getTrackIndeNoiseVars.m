function [trackNoiseVars] = getTrackIndeNoiseVars (kalmanInfoLink,tracksFinal,trackIdx,schemeIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;


trackNoiseVars=NaN(size(tracksFinal(trackIdx).tracksFeatIndxCG,2),4);
for t = startTime:endTime
    featureIdx=tracksFinal(trackIdx).tracksFeatIndxCG(tIdx);
    
    if featureIdx ~=0
        trackNoiseVars(tIdx,:)=diag(kalmanInfoLink(t).noiseVar(:,:,featureIdx,schemeIdx));
    else
        trackNoiseVars(tIdx,:)=  zeros(1,4);
    end
    tIdx=tIdx+1;
end
