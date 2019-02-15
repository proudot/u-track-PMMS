function tracksOut = reformatTracks(trackFeatInfo,detectInfo)
%SEPARATETRACKS reformats compound tracks into one structure per compound track
%
% tracksOut = separateTracks(trackFeatInfo)
% tracksOut = separateTracks(trackFeatInfo,detectionInfo)
%
%   Reformats the output of the U-Track package into separate structures,
%   one per track, with all associated info. Will also incorporate
%   additional info from any detection structure that is input, so long as
%   it follows the format of one structure per image, with each field
%   having a number of elements equal to the number of detected features in
%   that image. All fields will be copied over except those which have
%   identical names to the tracking fields.
%
%
%
% Hunter Elliott
% 8/2013
% Copied / modified from TrackingProcess.m

if nargin < 2
    detectInfo = [];
end


            
% Determine the number of compound tracks
nCompoundTracks = zeros(numel(trackFeatInfo), 1);
for i = 1:numel(trackFeatInfo)
    nCompoundTracks(i) =  size(trackFeatInfo(i).tracksCoordAmpCG,1);
end
nTracksTot = [0 cumsum(nCompoundTracks(:))'];

% Fail fast if no track
if sum(nCompoundTracks) == 0
    tracksOut = struct.empty(1,0);
    return        
end


%Track field names and corresponding column indices
trackFields = {'x','y','z','A','x_pstd','y_pstd','z_pstd','A_pstd'};
fInd   = [ 1   2   3    4     5      6      7       8    ];
nField = numel(trackFields);

if isempty(detectInfo)
    detFields = {};        
    nDetField = 0;
else
    
    %Get fields from input detection structure
    detFields = fieldnames(detectInfo);
    nDetField = numel(detFields);
    %Remove duplicates and mis-sized fields
    isBad = false(nDetField,1);
    for j= 1:nDetField
        isBad(j) = any(strcmp(detFields{j},trackFields)) || ...
            ~isequal(numel(detectInfo(1).(detFields{j})),numel(detectInfo(1).x));%Lazy way to check size by comparing with x field, assumed to be present

    end
    detFields = detFields(~isBad);
    nDetField = numel(detFields);    
end

strCell = vertcat([trackFields detFields'],arrayfun(@(x)([]),1:(nField+nDetField),'Unif',0));%Cell array for initializing structure
tracksOut(sum(nCompoundTracks),1) = struct(strCell{:});

for i = find(nCompoundTracks)'
    % Get the x and y coordinate of all compound tracks
    for  j = 1 : nCompoundTracks(i)
        iTrack = nTracksTot(i) + j ;
        tracksOut(iTrack).number = i;
        for k= 1:nField        
            tracksOut(iTrack).(trackFields{k}) = trackFeatInfo(i).tracksCoordAmpCG(j, fInd(k):8:end);                        
        end
        if isfield(trackFeatInfo, 'label'),
            tracksOut(iTrack).label = trackFeatInfo(i).label;
        end
    end

    % Fill split events NaNs in compound tracks
    nTimes = numel(tracksOut(iTrack).(trackFields{1}));
    splitEvents = find(trackFeatInfo(i).seqOfEvents(:,2)==1 & ~isnan(trackFeatInfo(i).seqOfEvents(:,4)))';
    eventTimes = trackFeatInfo(i).seqOfEvents(splitEvents,1)-1;
    for iEvent = splitEvents(eventTimes < nTimes)
        iTrack1 = nTracksTot(i)+ trackFeatInfo(i).seqOfEvents(iEvent,3);
        iTrack2 = nTracksTot(i)+ trackFeatInfo(i).seqOfEvents(iEvent,4);
        t = trackFeatInfo(i).seqOfEvents(iEvent,1)-1;
        for k = 1:nField
            tracksOut(iTrack1).(trackFields{k})(t) = tracksOut(iTrack2).(trackFields{k})(t);
        end                
    end

    % Fill merge events NaNs in compound tracks
    mergeEvents = find(trackFeatInfo(i).seqOfEvents(:,2)==2 & ~isnan(trackFeatInfo(i).seqOfEvents(:,4)))';
    eventTimes = trackFeatInfo(i).seqOfEvents(mergeEvents,1)-1;
    for iEvent = mergeEvents(eventTimes < nTimes)
        iTrack1 = nTracksTot(i)+ trackFeatInfo(i).seqOfEvents(iEvent,3);
        iTrack2 = nTracksTot(i)+ trackFeatInfo(i).seqOfEvents(iEvent,4);
        t = trackFeatInfo(i).seqOfEvents(iEvent,1);
        for k = 1:nField
            tracksOut(iTrack1).(trackFields{k})(t) = tracksOut(iTrack2).(trackFields{k})(t);
        end       
    end
    
    %Extract corresponding values from detection structure
    iFrameCurr = trackFeatInfo(i).seqOfEvents(1,1):trackFeatInfo(i).seqOfEvents(end,1);%Frames for current track;
    iPres = find(trackFeatInfo(i).tracksFeatIndxCG~=0);
    iFramePres = iFrameCurr(iPres);
    for k = 1:nDetField        
        tracksOut(iTrack).(detFields{k}) = nan(1,numel(iPres));    
        for j = 1:numel(iFramePres)
            tracksOut(iTrack).(detFields{k})(iPres(j)) = detectInfo(iFramePres(j)).(detFields{k})(trackFeatInfo(i).tracksFeatIndxCG(iPres(j)));                        
        end                
    end    
end




