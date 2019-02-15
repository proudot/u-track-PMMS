function [trackStateCovs] = getTrackStateCovs (kalmanInfoLink,tracksFinal,trackIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;

trackStateCovs=NaN(4,4,size(tracksFinal(trackIdx).tracksFeatIndxCG,2));
for t = startTime:endTime
    trackStateCovs(:,:,tIdx)=kalmanInfoLink(t).stateCov(:,:,tracksFinal(trackIdx).tracksFeatIndxCG(tIdx));
    tIdx=tIdx+1;
end