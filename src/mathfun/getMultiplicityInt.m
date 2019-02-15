% [r, udata, sdata] = getMultiplicityInt(data,index,skipFiniteCheck) returns the occurrences/Multiplicity of the elements of 'data'
%
% Inputs:
%         data : n-dimensional input array (Nx1 linearized)
%         index : (optional) 1) boolean, default: false
%                               use a unique index to minimize memory usage
%                            2) unique data values
%                            3) The third output of unique, ic
%         skipFiniteCheck : (optional) boolean, default: false
%                            Save processing when using double
%                            by skipping check for NaN and Inf
%                            
%
% Outputs: 
%          rep : # of occurrences for each element of 'data'
%        udata : sorted 1-D array of unique values in 'data'
%        sdata : sorted 1-D array of values in 'data'
%
% Note: NaN/Inf elements in input data are ignored
%
% This is an optimization for getMultiplicity for arrays of integers in a 
% narrow range. The optimization consumes more memory with a larger range of values
% in data and requires that the values be integers, or mappable to integers
%
% Essentially, this function expands the domain for counting using accumarray
% to be as large as possible: anything shiftable to between 1 and 2^64-1*.
%
% Significant time is saved if given an integer data type since the function does
% not have to check for non-finite values.
%
% The second input, if not false,  uses unique or the output of unique to minimize
% memory consumption at the cost of processing time.
%
%
% Alternatives:
% 1) getMultiplicity(data)
% 2) histc(data,unique(data))
% 3) accumarray(data(:),1)
%
% *MATLAB can only address arrays up to 2^48-1 as of 2014/09/13

% Mark Kittisopikul, 2014/08/11

function [rep, udata, sdata] = getMultiplicityInt(data,index,skipFiniteCheck)

data = data(:);
udata = [];
sdata = [];
dataTypeMax = Inf;
dataTypeMin = -Inf;
dataType = class(data);

if(isempty(data))
    rep = [];
    return;
end

if(nargin < 2)
    index = false;
end

if(nargin < 3)
    skipFiniteCheck = false;
end

if(isfloat(data))
    % only check for NaN or Inf if not an integer datatype
    if(~skipFiniteCheck)
        % Doing this consumes more time than the rest of the function
        data = data(isfinite(data));
    end
    % flintmax availabe as of MATLAB 8.1 (R2013a)
    % 2^53 for double, 2^24 for single
    dataTypeMax = flintmax(dataType);
    dataTypeMin = realmin(dataType);
elseif(isinteger(data))
    dataTypeMax = intmax(dataType);
    dataTypeMin = intmin(dataType);
elseif(islogical(data))
    dataTypeMax = 1;
    dataTypeMin = 0;
end

if(isscalar(index))
    if(index)
        [udata, ia, ic] = unique(data);
         if(nargout > 2)
            sdata = data(ia);
         end
        data = ic;
    end
elseif(length(index) == length(data))
    % if you have the third output of unique, then just pass that as data
    % index is ic
    if(nargout > 1)
        udata = data(index);
     end
     data = index;
else
    % index is udata
    assert(strcmp(class(index),dataType));
    udata = index;
    map = sparse(max(udata),1);
    map(udata) = 1:length(udata);
    data = map(data);
    data = full(data);
end


% set the minimum value to 1 to make the maximum value small
if(data(1) == dataTypeMin)
    minValue = dataTypeMin;
else
    minValue = min(data);
end
if(minValue <= 0 || nargout > 1)
    maxValue = max(data);
    if(minValue <= 0 && maxValue == dataTypeMax || maxValue-minValue >= dataTypeMax)
        if(strcmp(dataType,'uint64') || strcmp(dataType,'int64'))
            err = MException('getMultiplicityInt:cannotShiftMinToOne', ...
                ['Cannot shift minimum to one: ', ...
                'Argument contains %.0f, ', ...
                'and the maximum value, %.0f, for type %s.' ], ...
                minValue, maxValue, class(data));
            throw(err);
        elseif(minValue < 0)
            % Could map into uint64, but 63 bits should be good enough
            data = int64(data);
        else
            data = uint64(data);
        end
    end
end
% avoid add/substract by zero for efficiency
if(minValue ~= 1)
    data = data - double(minValue) + 1;
end

try
    %where the magic happens
    rep = accumarray(data,1);
catch err
    if(strcmp(err.identifier,'MATLAB:accumarray:nonPosIntIndValues'))
        err2 = MException('getMultiplicityInt:nonIntegerValues', ...
            ['Argument, data, has non-integer values. ',  ...
             'Use getMultiplicity instead.']);
        err = addCause(err,err2);
        rethrow(err);
    elseif(strcmp(err.identifier,'MATLAB:accumarray:pmaxsize') || ...
           strcmp(err.identifier,'MATLAB:nomem'))
        % maxValue - minValue is too big
        if(isscalar(index) && ~index)
            % attempt to index
            warning(['getMultiplicityInt: Range of values in data is large. ' ...
                     'Set second parameter to true or use getMultiplicity ' ...
                     'for faster execution.'])
            out = cell(1,3);
            [out{1:max(nargout,1)}] =  getMultiplicityInt(data + minValue - 1,true,skipFiniteCheck);
            [rep,udata,sdata] = out{:};
            return;
        else
            % if we already tried to index, then fail
            err2 = MException('getMultiplicityInt:largeRangeofValues', ...
                ['The range of values in data is too large. ', ...
                 'Use getMultiplicity instead or index using the second parameter.']);
            err = addCause(err,err2);
            rethrow(err);
        end
    else
        rethrow(err);
    end
end


% if unique values requested, then return those with nonzero repeitions
if(nargout > 1 && isempty(udata))
    udata = minValue : maxValue;
    udata = udata(rep ~= 0);
end

if(nargout > 2 && isempty(sdata))
    % do run length decompression like rude in File Exchange
    repcum = cumsum(rep(rep ~= 0));
    idx = data;
    idx(1) = 1;
    idx(2:end) = 0;
    idx( repcum(1:end-1) + 1) = 1;
    idx = cumsum(idx);
    sdata = udata(:,idx);
end

rep = rep(rep ~= 0)';

end
