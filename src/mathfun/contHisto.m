function varargout = contHisto(data,distribution,normalized,verb,n)
%CONTHISTO plots data in a continuous histogram (good for plotting data with error margins)
%
%SYNOPSIS [Y,X,YDetail] = contHisto(data,distribution,normalized,verb,n)
%
%INPUT    data        : vector containing the values x+-s in pairs 
%                       of the form [x s]
%         dist        : (opt) {'norm'} or 'unif' depending on whether the
%                       data should be represented by normal or uniform 
%                       distributions
%         normalized  : (opt) {1}/0 if the surface each data point
%                             contributes should be normalized
%         n           : (opt) number of points along the x-axis (500). min 10
%         verb        : (opt) {0}/1 whether the histogram should be plotted.
%
%         If the function is called with no output arguments, the histogram
%         is plotted, anyway. Norm-plots are plotted using the 'area'
%         function, ident-plots with 'stairs'
%
%OUPTUT   Y       : vector of length n ('norm') or 2xlength data ('unif') with
%                   the y-data for the plot 
%         X       : vector with the corresponding x-data
%         YDetail : n-by length(x) array with data for each entry
%                   individually
%
%c: 12/03, jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign defaults
plottype      = 'norm'; %normOrIdent
normalizeData = 1; %normalized
verbose       = 0; %verb
xLength       = 500; %n

%------TEST INPUT------

%data
if nargin == 0 || isempty(data)
    error('not enough input arguments for contHisto')
elseif size(data,2) ~= 2 
    error('data must be a N-by-2 vector')
end

%distribution
if nargin < 2 || isempty(distribution)
    %we take the default
elseif strcmp(distribution,'norm') || strcmp(distribution,'unif')
    plottype = distribution;
else
    error('bad option for distribution');
end

%normalized
if nargin < 3 || isempty(normalized)
    %we take the default
elseif all(normalized == 1) || all(normalized == 0)
    normalizeData = normalized(1);
else
    error('bad option for normalized');
end

%verbose
if nargin < 4 || isempty(verb)
    %we take the default
elseif all(verb == 1) || all(verb == 0)
    verbose = verb(1);
else
    error('bad option for verb');
end

%n
if nargin < 5 || isempty(n)
    %we take the default
elseif n>10
    xLength = n;
else
    error('n has to be a positive integer > 10')
end


%------END TEST INPUT

%------TEST OUTPUT
if nargout == 0
    verbose = 1;
end
%--END TEST OUTPUT


%-----CALCULATE PLOT DATA

%switch according to plottype

switch plottype
    
    case 'unif'
        %build y out of summed uniformly distributed blocks
        
        %1) get interval borders, weight
        intervalStart = diff(data,1,2); %X-S
        intervalEnd   = sum(data,2); %X+S
        nData = size(data,1);
        
        if normalizeData
            weights = 1./(2*data(:,2));
        else
            weights = ones(nData,1);
        end
        
        %2) generate lists of xData with the corresponding entry numbers.
        sList = [[1:nData]',intervalStart];
        eList = [-[1:nData]',intervalEnd];
        
        xList = sortrows([sList;eList],2);
        
        %3) assign Y-Data.
        YDetail = zeros(2*nData,nData);
        
        for i = 1:nData
            yIndex = [find(xList(:,1)==i):find(xList(:,1)==-i)]; %find the range spanned by the current dataPoint
            YDetail(yIndex,i) = weights(i);
        end
        
        %calc X and Y
        X = xList(:,2);
        Y = sum(YDetail,2);
        
        %plot if selected
        if verbose
            figure;
            stairs(X,Y);
        end
        
        %end case 'ident'
        
    case 'norm'
        %build y out of individual normal pdfs
        
        %1) calculate X, mu and sigma - arrays; assign weights
        xMin = min(data(:,1)-2*data(:,2));
        xMax = max(data(:,1)+2*data(:,2));
        
        X = [xMin:((xMax-xMin)/(xLength-1)):xMax]';
        nData = size(data,1);
        
        XMat = repmat(X,[1,nData]);
        
        muMat = repmat(data(:,1)',[xLength,1]);
        siMat = repmat(abs(data(:,2)'),[xLength,1]);
        
        %to sort of have an idea of how many counts we have, normalize the
        %data with the mean hight of each distribution
        invWeights = sqrt(pi*2)*data(:,2)';
        standardHeightMultiplicator = 1/mean(1./invWeights);
        
        if normalizeData %normal distribution is normalized already!
            weights = ones(xLength,nData)*standardHeightMultiplicator;
        else %we use the weights to un-normalize the data
            weights = repmat(invWeights,[xLength,1]); %max is 1 for all, no multiplications!
        end
        
        %calculate y and yDetails
        YDetail = normpdf(XMat,muMat,siMat).*weights;
        
        Y = sum(YDetail,2);
        
        %plot if selected
        if verbose
            figure;
            area(X,Y);
            figure
            area(X,YDetail);
        end
end

%assign varargout
if nargout > 0
    varargout{1} = Y;
end
if nargout > 1
    varargout{2} = X;
end
if nargout > 2
    varargout{3} = YDetail;
end

        
        