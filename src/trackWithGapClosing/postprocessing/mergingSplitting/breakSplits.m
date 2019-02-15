function tracksFinal = breakSplits(tracksFinal,split2remove,iTrack)
%BREAKSPLITS breaks selected splits in a tracks-structure
%
% SYNOPSIS: tracksFinal = breakSplits(tracksFinal,split2remove,iTrack)
%
% INPUT tracksFinal: tracksFinal structure from trackCloseGapsKalman.
%		split2remove: 2-by-n numeric array with
%           [motherSegment;daughterSegment] of splits to remove, or cell
%           array of length nTracks containing split2remove arrays (empty if
%           there is nothing to remove).      
%		iTrack : if split2remove is numeric, indicate the track number it
%           corresponds to 
%
% OUTPUT tracksFinal: tracksFinal structure where selected track has been
%           split. Its length is nOld+k, where k is the number of broken
%           splits  
%
% REMARKS breakSplits has only been tested on merge-less tracks, but it
%         may work on tracks including merges, too.
%
% created with MATLAB ver.: 7.8.0.8205 (R2009a) Beta (Mac Intel 64-bit) on Mac OS X  Version: 10.5.6 Build: 9G2141 
%
% created by: jonas
% DATE: 29-Apr-2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
if nargin < 2 || isempty(tracksFinal) || ~isstruct(tracksFinal)
    error('please call breakSplits with at least two input arguments and nonempty track-structure')
end

nTracks = length(tracksFinal);

% check split2remove
if isempty(split2remove)
    return
end
if isnumeric(split2remove) 
    if nargin < 3 || isempty(iTrack)
        error('please specify which track split2remove corresponds to by supplying the argument iTrack')
    end
    todoList = iTrack;
    split2remove = {split2remove};
elseif iscell(split2remove)
    todoList = find(~cellfun(@isempty,split2remove));
    split2remove = split2remove(todoList);
else
    error('split2remove has either to be a cell array or a numeric array')
end


%% main loop

for tt = 1:length(todoList)
    % this was originally part of another function
    iTrack = todoList(tt);
    div2remove = split2remove{tt};

% need to sort div2remove with descending d2r(2), so that there
        % won't be trouble with indices after deleting tracks.
        if size(div2remove,2) > 1
            div2remove = sortrows(div2remove',-2)';
        end
        for d2r = div2remove % d2r: iSeg,otherSegment
            
            
            % find other tracks that link to this track
            transferIdx = d2r(2);
            seg2check = d2r(2);
            while ~isempty(seg2check)
                % look for divisions off d2r(2)
                moreDivIdx = find(tracksFinal(iTrack).seqOfEvents(:,2) == 1 & tracksFinal(iTrack).seqOfEvents(:,4) == seg2check(1)); 
                seg2check(1) = [];
                if ~isempty(moreDivIdx)
                    % remember additional segment indices
                    transferIdx = [transferIdx;tracksFinal(iTrack).seqOfEvents(moreDivIdx,3)];
                    seg2check = [seg2check;tracksFinal(iTrack).seqOfEvents(moreDivIdx,3)];
                end
            end
            
            % transfer all 
            nTracks = nTracks + 1;
            tracksFinal(nTracks).tracksCoordAmpCG = tracksFinal(iTrack).tracksCoordAmpCG(transferIdx,:);
            tracksFinal(nTracks).tracksFeatIndxCG = tracksFinal(iTrack).tracksFeatIndxCG(transferIdx,:);
            soeTransfer = ismember(tracksFinal(iTrack).seqOfEvents(:,3),transferIdx);
            tracksFinal(nTracks).seqOfEvents = tracksFinal(iTrack).seqOfEvents(soeTransfer,:);
            % turn start of segment d2r(2) into a birth
            rmIdx = all(bsxfun(@eq,tracksFinal(nTracks).seqOfEvents(:,2:4),[1,d2r(2),d2r(1)]),2);
            tracksFinal(nTracks).seqOfEvents(rmIdx,4) = NaN;
            % update indices in new tracks
            for ii=1:length(transferIdx) % the ith entry in transferIdx becomes i
                iiIdx = find(tracksFinal(nTracks).seqOfEvents(:,3:4) == transferIdx(ii));
                iiIdx = iiIdx + 2*size(tracksFinal(nTracks).seqOfEvents,1); % shift index to last two cols
                tracksFinal(nTracks).seqOfEvents(iiIdx) = ii;
            end
            % remove starting NaNs (or zeros) in new tracks
            firstGood = find(any(tracksFinal(nTracks).tracksFeatIndxCG>0,1),1,'first');
            tracksFinal(nTracks).tracksFeatIndxCG(:,1:firstGood-1) = [];
            tracksFinal(nTracks).tracksCoordAmpCG(:,1:8*(firstGood-1)) = [];
            % remove trailing NaNs in new tracks. Update track-end
            lastGood = find(any(tracksFinal(nTracks).tracksFeatIndxCG>0,1),1,'last');
            tracksFinal(nTracks).tracksFeatIndxCG(:,lastGood+1:end) = [];
            tracksFinal(nTracks).tracksCoordAmpCG(:,8*(lastGood)+1:end) = [];
            tracksFinal(nTracks).seqOfEvents(:,1) = min(tracksFinal(nTracks).seqOfEvents(:,1),...
                min(tracksFinal(nTracks).seqOfEvents(:,1)) + lastGood - 1);
          
            % remove track in old data
            tracksFinal(iTrack).tracksCoordAmpCG(transferIdx,:) = [];
            tracksFinal(iTrack).tracksFeatIndxCG(transferIdx,:) = [];
            tracksFinal(iTrack).seqOfEvents(soeTransfer,:) = [];
           % remove trailing NaNs in oldTracks
                lastGood = find(any(tracksFinal(iTrack).tracksFeatIndxCG>0,1),1,'last');
            tracksFinal(iTrack).tracksFeatIndxCG(:,lastGood+1:end) = [];
            tracksFinal(iTrack).tracksCoordAmpCG(:,8*(lastGood)+1:end) = [];
            tracksFinal(iTrack).seqOfEvents(:,1) = min(tracksFinal(iTrack).seqOfEvents(:,1),...
                min(tracksFinal(iTrack).seqOfEvents(:,1)) + lastGood - 1);
             % update indices in old data  
            for idx = unique(tracksFinal(iTrack).seqOfEvents(:,3))'
                delta = sum(idx > transferIdx);
                if delta > 0
                    iiIdx = find(tracksFinal(iTrack).seqOfEvents(:,3:4) == idx);
                    % shift like above
                    iiIdx = iiIdx + 2 * size(tracksFinal(iTrack).seqOfEvents,1);
                    tracksFinal(iTrack).seqOfEvents(iiIdx) = idx - delta;
                end
            end
        end
end