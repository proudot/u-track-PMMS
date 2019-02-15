function randStat = randStatisticCorrected(groups1,groups2)
%RANDSTATISTICCORRECTED calculates the corrected rand for two grouping procedures
%
% SYNOPSIS: randStat = randStatisticCorrected(groups1,groups2)
%
% INPUT groups: vectors containing the group labels according to procedure
%               1 and 2, respectively, for every data item. Group labels
%               have to be finite integers > 0
%
% OUTPUT randStat statistic. 1 if perfect agreement, 0 if random grouping,
%        negative if actively messed up 
%
% REMARKS see e.g. Jain and Dubes (1984) Algorithms for Clustering Data.
%         Prentice Hall, New Jersey. pp 172-177
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 13-Sep-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%==================
%% TEST INPUT
%==================

% only check for the nubmer of input arguments - all other tests will be
% performed by frequencyTable.m
if nargin < 2
    error('randStatisticCorrected needs two nonempty input arguments!');
end

%==================

%============================
%% GENERATE FREQUENCY TABLE
%============================

% try/catch because we want to rename the error
try
    table = frequencyTable(groups1,groups2);
catch
    % if any error: read message, replace "frequencyTable" with our
    % functionname, and rethrow the error
    err = lasterror;
    
    % check for the proper ID, and change message
    if strmatch('FREQUENCYTABLES:WRONGINPUT',err.identifier)
       err.message = regexprep(err.message,'frequencyTable',mfilename);
    end
    
    % rethrow the error
    rethrow(err)
end

%===============================


%===============================
%% CALCULATE CORRECTED RAND
%===============================

% calculate the binomial coefficients. Since they are all (n 2), we just
% take n*(n-1)/2

% sum over all table frequences
bsumNij = sum(table(:).*(table(:)-1))/2;

% sum over the row frequences
rowSum = sum(table,2);
bsumNi = sum(rowSum.*(rowSum-1))/2;

% sum over col frequences
colSum = sum(table,1);
bsumNj = sum(colSum.*(colSum-1))/2;

% number of possible pairs
nData = length(groups1);
nPairs = nData *(nData-1)/2;

% calculate corrected rand
randStat = (bsumNij - bsumNi*bsumNj/nPairs)/((bsumNi+bsumNj)/2 - bsumNi*bsumNj/nPairs);

