function [trackStateCovs] = getTrackIndeStateCovs (kalmanInfoLink,tracksFinal,trackIdx,schemeIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;


trackStateCovs=NaN(4,4,size(tracksFinal(trackIdx).tracksFeatIndxCG,2));
for t = startTime:endTime
    featureIdx=tracksFinal(trackIdx).tracksFeatIndxCG(tIdx);
    
    if featureIdx ~=0
        trackStateCovs(:,:,tIdx)=kalmanInfoLink(t).stateCov(:,:,featureIdx,schemeIdx);
    else
        trackStateCovs(:,:,tIdx)=  zeros(4);
    end
    tIdx=tIdx+1;
end
