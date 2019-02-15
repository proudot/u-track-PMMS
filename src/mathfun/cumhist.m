function [no,xo] = cumhist(y,x)
%CUMHIST  Histogram.
%   N = CUMHIST(Y) bins the elements of Y into 10 equally spaced containers
%   and returns the sum of the number of elements in the first i containers (i=1..10).
%   If Y is a matrix, HIST works down the columns.
%
%   N = CUMHIST(Y,M), where M is a scalar, uses M bins.
%
%   N = CUMHIST(Y,X), where X is a vector, returns the distribution of Y
%   among bins with centers specified by X. 
%
%   [N,X] = CUMHIST(...) also returns the position of the bin centers in X.
%
%   CUMHIST(...) without output arguments produces a histogram bar plot of
%   the results.
%
%   See also HIST, HISTC.

% c: 31/5/01	dT

if nargin == 0
    error('Requires one or two input arguments.')
end
if nargin == 1
    x = 10;
end

if min(size(y))==1, y = y(:); end
if ischar(x) || ischar(y)
    error('Input arguments must be numeric.')
end

if isempty(y),
    if length(x) == 1,
       x = 1:x;
    end
    n = zeros(size(x)); % No elements to count
else
    if length(x) == 1
        miny = min(min(y));
        maxy = max(max(y));
    	  if miny == maxy,
    		  miny = miny - floor(x/2) - 0.5; 
    		  maxy = maxy + ceil(x/2) - 0.5;
     	  end
        binwidth = (maxy - miny) ./ (x-1);
        x = miny + binwidth*(0:(x-1));
    end
  
    for i=1:length(x)
    n(i)=nnz(find(y<=x(i)));
    end;
end

if nargout == 0
    plot(x,n);
else
    xo=x;
    no=n;
end
