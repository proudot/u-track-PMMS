function cutValue = splitModes(data,loHi,jumpMax,verbose)
%SPLITMODES cut histogram where there is no or little data
%
% SYNOPSIS cutValue = splitModes(data,loHi,jumpMax)
%
% INPUT    data : data to split
%          loHi : (opt) number of minima that should be checked to the left
%                 and right of the initial guess for the cutoff. The lowest
%                 minimum right of the highest peak will be chosen.
%                 Default: [1,2] - 1 to the left, 2 to the right. Note that
%                 the number of minima does not include the minumum closest
%                 to the initial guess, i.e. [0,0] chooses the minimum
%                 closest to the initial guess without looking for more.
%          jumpMax : (opt) if true, minimum is searched to the left of the
%                 highest maximum, too. Default: 0
%          verbose : (opt) if 1, output will be plotted. Default: 0
%
% c: 07/2008 Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(loHi)
    loHi = [1,2];
end
if nargin < 3 || isempty(jumpMax)
    jumpMax = false;
end
if nargin < 4 || isempty(verbose)
    verbose = false;
end

% first guess via cutFirstHistMode - put as 1 to see histograms
if iscell(data)
    % data is already in counts/bins
    [cutIdx,cutVal] = cutFirstHistMode(data{:},0);
    % create the spline
    sp = spline(data{2},data{1});
else
[cutIdx, cutVal,sp] = cutFirstHistMode(data,0);
end

% now check the local minima in the vicinity of the cutoff
spder = fnder(sp);
zeroList = fnzeros(spder);
zeroList = zeroList(1,:);
% evaluate
zeroVals = fnval(sp,zeroList);

% look in zeroList. Find one value before cutVal, three after. Go into
% zeroVals and find lowest minimum

%KJ: this original line assumes that the index closest to the first guess
%is a minimum, but that is not necessarily the case
%so I am modifying it to make sure that what is designated as the closest
%index is for sure a minimum
%this is particularly important when hiLo = [0 0]
%
% [dummy,closestIdx] = min(abs(zeroList - cutVal));

distFromCutVal = zeroList - cutVal;
closestIdxNeg = max(length(distFromCutVal(distFromCutVal<0)),1);
closestIdxPos = closestIdxNeg + 1;
closestIdxBoth = [closestIdxNeg closestIdxPos];
zeroValsBoth = zeroVals(closestIdxBoth);
[dummy,minIdx] = min(zeroValsBoth);
closestIdx = closestIdxBoth(minIdx);

% check only the minima that are close by; two to the right and
% one to the left (don't forget that between two minima there will
% always be a maximum!)
indexList = (closestIdx-loHi(1)*2):(closestIdx + loHi(2)*2);
indexList(indexList < 1 | indexList > length(zeroVals)) = [];

% also, never allow going left of the highest maximum!
[maxVal,maxValIdx] = max(zeroVals);
if ~jumpMax
    indexList(indexList < maxValIdx) = [];
end

% find lowest
[dummy, cutIdx] = min(zeroVals(indexList));
% and determine break value
cutValue = zeroList(indexList(cutIdx));

% if verbose: plot.
if verbose
    figure
    ah = gca;
    if iscell(data)
        bar(data{2},data{1})
    else
    optimalHistogram(ah,data,1,0);
    end
    hold on
    plot(ah,[cutVal;cutVal],[0,maxVal],':r')
    plot(ah,[cutValue;cutValue],[0,maxVal],'r')
end
