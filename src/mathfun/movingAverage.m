function [outputVector, result] = movingAverage (varargin)
% movingAverage calculates a moving average over the input value specified
% by method 'method'
%
% SYNOPSIS       [outputVector, result] = movingAverage (inputVector, windowSize, method)
%
% INPUT          inputVector : the vector with input data (1,1:max)
%                windowSize : the window size used to average (default 5)
%                method : the method used to average (default median)
%                
% OUTPUT         outputVector : vector with averaged data
%                result : result of the operation: 0 = ok, 1 = error
%
% DEPENDENCIES   movingAverage   uses {  }
%                                  
%                movingAverage  is used by {  }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Sep 04          Initial release

if nargin < 1
    fprintf (1, 'Error: movingAverage: not enough input parameters!\n');
end

% Get input data
l=length(varargin);
switch l
    case 1
        inputVector = varargin{1};
        windowSize = 5;
        method = 'median';
    case 2
        inputVector = varargin{1};
        windowSize = varargin{2};
        method = 'median';
    case 3
        inputVector = varargin{1};
        windowSize = varargin{2};
        method = varargin{3};
    otherwise
        fprintf (1, 'Error: movingAverage: Too many input parameters!\n');
        return;
end

% Determine maximum loop count
maxLoopCount = length(inputVector) - windowSize + 1;

% Initialize output data matrix
workMatrix = zeros (windowSize, maxLoopCount);

for iCount = 1 : maxLoopCount
    % Store in the work matrix
    matrixEntry = inputVector(iCount:iCount+windowSize-1);
    workMatrix(:,iCount) = matrixEntry';
end

% Depending on the method do the appropriate filtering
switch lower(method)
   case 'median'
      outputVector = medfilt1 (inputVector, windowSize);
      result = 0;
   otherwise
      fprintf (1, ['Method ' method ' not implemented yet.']);
      outputVector = [];
      result = 1;
end