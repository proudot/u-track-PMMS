function [trackBestSchemes] = getTrackIndeBestSchemes (kalmanInfoLink,tracksFinal,trackIdx)
% Right now algorithm only focus on the one tracklet tracks

sOE=tracksFinal(trackIdx).seqOfEvents;
startTime=sOE(sOE(:,2)==1,1);
endTime=sOE(sOE(:,2)==2,1);
tIdx=1;


trackBestSchemes=NaN(size(tracksFinal(trackIdx).tracksFeatIndxCG,2),1);
for t = startTime:endTime
    featureIdx=tracksFinal(trackIdx).tracksFeatIndxCG(tIdx);    
    if featureIdx ~=0
        trackBestSchemes(tIdx)=kalmanInfoLink(t).scheme(featureIdx,2);
    else
        trackBestSchemes(tIdx)=  0;
    end
    
    tIdx=tIdx+1;
end
