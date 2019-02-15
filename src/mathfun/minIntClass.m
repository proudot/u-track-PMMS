function minClass = minIntClass(n)
%MININTCLASS returns minimum integer class required to store specified number
%
% minClass = minIntClass(n);
%
%   This function returns the string name of the smallest bit-depth
%   unsigned integer class that can contain a number as high as the input
%   value n. If n <=1 it returns 'logical'
%
%Hunter Elliott
%6-26-2013

if n < 0
    error('This function doesn''t yet support signed integers! Feel free to add support for these!!')
end

nPow = log2(double(n)+1);

    
if  nPow <=1
    minClass = 'logical';
    return
elseif nPow > 64    
    error('n too large to be contained in integer class!')
end

bitPos = 8*2.^(0:3);

minClass = ['uint' num2str(bitPos(find(bitPos - nPow>=0,1,'first')))];

