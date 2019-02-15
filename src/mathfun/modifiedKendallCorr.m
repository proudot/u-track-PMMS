function  [tau,tauAmp] = modifiedKendallCorr(x,varargin)
%Modified Kendall cross-correlation coefficient
%Temporal information added 
%Amplitude information added
%TS with NaN
%Under construction (missing p-value)
%Based on "A New Measure of Rank Correlation", M.G.Kendall. Biometrika, Vol. 30, No 1/2, pp 81-93
%Marco Vilela, 2012

 

 

 
ip = inputParser;
ip.addRequired('x',@(x) isvector(x));
ip.addOptional('y',x,@(x) isvector(x));
ip.addParamValue('local',numel(x)-1,@isscalar);
ip.addParamValue('alpha',0.05,@isscalar);
ip.addParamValue('maxLag',0,@isscalar);

 

 
ip.parse(x,varargin{:});
local    = ip.Results.local;
alpha    = ip.Results.alpha;
y        = ip.Results.y;
maxLag   = ip.Results.maxLag;

 
x = x(:);
y = y(:);

 
if ~all( size(x) == size(y) )
    error('x and y have different size')
else
    nObs   = numel(x) - max( [sum(isnan(x)) sum(isnan(y))] );
    localN = nObs - 1;
    if nObs <= 6 %This number comes from a lookup table (Original Kendall Paper) - Hard(impossible) to assume normalitity for smaller sample size
        error('Insufficient number of points')
    end
end

 
if maxLag > 0 && (nObs - local) < maxLag
    disp('maximum lag is too large. Reducing it to the max allowed')
   % maxLag = max((nObs - local) - 1,0);
end

 

 
normalization       = local*(2*nObs - local - 1)/2;
contM               = getCorrMatrix(ones(size(x)),local,[]);
[contX,ampM(:,:,1)] = getCorrMatrix(x,local,contM);
[contY,ampM(:,:,2)] = getCorrMatrix(y,local,contM);

 
%Matrix with [-1 0 1] for mismatch, no-match and match of slopes
match = contX.*contY;
tau   = sum(match(:))/normalization;

 
%Matrix with mean of [delta(i,j)] 
meanAmp  = mean(abs(ampM),3);
matchAmp = match.*meanAmp;
tauAmp   = sum(matchAmp(:))/sum(meanAmp(:));

 

 
contYr        = contY;
contXl        = contX;
ampM(:,:,3:4) = ampM(:,:,1:2);% 1 is for x and 2 for y
tauR          = NaN(1,maxLag-1);
tauL          = NaN(1,maxLag-1);
tauAmpR       = NaN(1,maxLag-1);
tauAmpL       = NaN(1,maxLag-1);

 

 
for iLag = 1:maxLag
    %Change this code - create a template that change with lag
    currNorm = normalization - local*iLag;

    
    %Shifting Y to the past
    [contYr,ampM(:,:,[2 4]),tauR(iLag),tauAmpR(iLag)]...
             = calculateCorr(contYr,contX,ampM(:,:,[2 4]),currNorm);

    
    %Shifting X to the past
    [contXl,ampM(:,:,[1 3]),tauL(iLag),tauAmpL(iLag)]...
             = calculateCorr(contXl,contY,ampM(:,:,[1 3]),currNorm);

    
end

 
tau    = [fliplr(tauL) tau tauR];
tauAmp = [fliplr(tauAmpL) tauAmp tauAmpR];

 
end%END OF MAIN FUNCTION

 
function [contTS,ampTSdiff] = getCorrMatrix(TS,neigh,contM)

 

 
hTS = hankel(TS);hTS(1,:) = [];
hTS(:,neigh+1:end)        = [];

 
if isempty(contM)   
    contTS    = hTS;
else
    ampTSdiff = (hTS - repmat(TS(1:end-1),1,neigh)).*contM;
    %Eliminating NaN - NaN points will not interfere with the counting, but they will still be there for the lag shifting
    ampTSdiff(isnan(ampTSdiff)) = 0;
    %
    contTS    = sign(ampTSdiff);
end

 
end%END OF getCorrMatrix

 
function [contZ,ampM,tau,tauAmp] = calculateCorr(contZ,contA,ampM,currNorm)
%contZ - comes from the time series being shifted
%contA - reference time series

 
    %Updating counter matrix - shifting up and adding zeros to the last row
    contZ  = circshift(contZ,-1);contZ(end,:) = 0;

 
    %Updating amplitude matrix 
    ampM(:,:,2) = circshift(ampM(:,:,2),-1);ampM(end,:,2) = 0;

    
    %Finding matches
    match    = contA.*contZ;

    
    %Calculating the mean amplitude of the slopes
    meanAmp  = mean( abs( ampM ),3);
    %Setting the sign for each amplitude slope
    matchAmp = match.*meanAmp;
    %
    tauAmp   = sum(matchAmp(:))/sum(meanAmp(:));
    tau      = sum(match(:))/currNorm;

    
end