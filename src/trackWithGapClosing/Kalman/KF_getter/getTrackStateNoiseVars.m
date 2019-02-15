function [trackStateNoiseVars] = getTrackStateNoiseVars (kalmanInfoLink,tracksFinal,trackIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;

trackStateNoiseVars=NaN(4,4,size(tracksFinal(trackIdx).tracksFeatIndxCG,2));
for t = startTime:endTime
    trackStateNoiseVars(:,:,tIdx)=kalmanInfoLink(t).noiseVar(:,:,tracksFinal(trackIdx).tracksFeatIndxCG(tIdx));
    tIdx=tIdx+1;
end