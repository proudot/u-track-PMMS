function table = frequencyTable(groups1, groups2)
%FREQUENCYTABLE generates frequency tables from two lists of labelings
%
% SYNOPSIS: table = frequencyTable(groups1, groups2)
%
% INPUT  groups: vectors with group number for every data item. Group
%                numbers should be finite integers > 0.
%
% OUTPUT table: matrix, where element m(i,j) contains the number of data
%               items that fall both into group i from groups1, and group j
%               from groups2  
%
% EXAMPLE  frequency of remainders of division by 3 and 5 among the even
%          numbers between 1 and 100
%
%           x=2:2:100;
%           m5=mod(x,5)+1; % add 1 because m5 would otherwise be 0-4
%           m3=mod(x,3)+1; % add 1, too
%           table = frequencyTable(m3,m5)
%           table =
%                3     4     3     3     3
%                4     3     3     3     4
%                3     3     4     4     3
%
%           The labels are sorted, so table(1,1) is the number of even
%           integers between 1 and 100 that are divisible by 3 and 5
%           without remainder.
%           
%           sum(table,1) and sum(table,2) give the row and column sums, and
%           sum(table(:)) == length(x) == length(m5) == length(m3)
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 13-Sep-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%======================
%% TEST INPUT
%======================

% use error ID, so that we can reuse the error messages in
% randStatisticCorrected

% all we want are two nonempty, finite vectors of the same length
if nargin < 2 || isempty(groups1) || isempty(groups2)
    error('FREQUENCYTABLES:WRONGINPUT',...
        'frequencyTables needs two nonempty input arguments')
end

% make vectors
groups1 = groups1(:);
groups2 = groups2(:);

% check length
if length(groups1) ~= length(groups2)
    error('FREQUENCYTABLES:WRONGINPUT',...
        'there must be the same number of elements in groups1 and groups2')
end

% combine to list of subscripts
subList = [groups1, groups2];

% check for finiteness
if ~all(all(isfinite(subList),1),2)
    error('FREQUENCYTABLES:WRONGINPUT',...
        'group labels cannot be inf or NaN!')
end

if any(any(subList < 1 | round(subList) ~= subList,1),2)
    error('FREQUENCYTABLES:WRONGINPUT',...
        'group labels have to be positive integers!')
end


%======================


%===========================
%% GENERATE FREQUENCY TABLE
%===========================

% every group index corresponds to a row- or col-idx. Thus, we can consider
% [groups1(i), groups2(j)] as subscripts into the output table. 
% Count the number of subscripts, convert the subscripts into a linera
% index and assign counts.

% preassign output
sizeTable = max(subList,[],1);
table = zeros(sizeTable);

% count entries in table
[sortedSubscripts, counts] = countEntries(subList,1);

% convert to linear index
linearIdx = sub2ind(sizeTable,sortedSubscripts(:,1),sortedSubscripts(:,2));

% assign table
table(linearIdx) = counts;