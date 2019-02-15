function pValue = chiSquareTest(contingencyTable)
%CHISQUARETEST tests whether frequencies observed in experiments are similar
%
% SYNOPSIS pValue = chiSquareTest(contingencyTable)
%
% INPUT    contingencyTable: r-by-c array of observed frequencies.
%               every row stands for an independent experimental condition
%                  example: r1: wildtype, r2: mutant1, r3: mutant2
%               every column stands for mutually exclusive observed
%               frequencies.
%                  example: c1: fast, c2: slow, c3: stationary
%
% OUTPUT   pValue: r-by-r matrix of probabilities.
%               The first column will be the comparison of all the rows
%               with the first condition, the second column will be the
%               comparison of all the rows with the second condition etc.
%               The p-values are given for the null hypothesis that the
%               observed frequencies for the two conditions are the same.
%               (the matrix is symmetric!)
%
% REMARKS  For 2-by-2 tables, you need at least a total of 20 observation.
%               (The test cannot be trusted otherwise!)
%               If the number of observations is between 20 and 40, no
%               frequency should be lower than 5. If the number of
%               observations is above 40, no frequency should be lower than
%               1. The code currently does not test for this - it is
%               assumed that the user takes care of this. 
%          The testing is done according to test 16a/b in "Parametric and
%          Nonparametric Statistical Procedures" by D.J.Sheskin. Method 2
%          of 16.8 is used.
%
%
% POSSIBLE EXTENSION  Use Fisher's exact test for frequencies lower than 20
%                     Write power of test in other half of matrix
%           
% c: 5/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=================
% TEST INPUT
%=================

% only check nargin
if nargin == 0 || isempty(contingencyTable)
    error('chiSquareTest needs a non-empty input!')
end

% get and check size of matrix
[rows, cols] = size(contingencyTable);
if any([rows,cols]==1)
    error('need at least a 2-by-2 contingency table for comparisons!')
end

%==================

%==========================
% LOOP TO MAKE COMPARISONS
%==========================

% First, we calculate the overall expected frequencies for all the tests.
% In general, the expected frequency of element (r,c) is calculated as the
% sum of all observations in row r multiplied by the sum of all
% observations in column c, divided by the total number of observations.

overallNumObservations = sum(contingencyTable(:));
overallExpectedFrequencies = ...
    repmat(sum(contingencyTable,2),1,cols) .*...
    repmat(sum(contingencyTable,1),rows,1) ./...
    overallNumObservations;

% Now loop through every pair of rows to compare their frequencies. The
% resulting pValue-matrix will be symmetric, as comparing a to b gives the
% same result as comparing b to a.

% assign matrix of chi2 values that will then be transformed into p-values
chi2Values = zeros(rows, rows);
for ri = 1:rows - 1
    for rj = ri + 1:rows
        % make table of the two rows to compare
        miniTable = contingencyTable([ri,rj],:);
        miniNumObservations = sum(miniTable(:));
        % only 2 rows for expected frequencies of mini table
        miniExpectedFrequencies = ...
            repmat(sum(miniTable,2),1,cols) .*...
            repmat(sum(miniTable,1),2,1) ./...
            miniNumObservations;
        
        % chiSquare is the sum of (Oij-Eij)^2/E'ij, where Oij is the observed
        % frequency of element ij, E is the expected frequency for the
        % mini-table, and E' is the expected frequency of the overall table
        % for that same element (with slightly different indices, of course)
        
        residualTable = ...
            (miniTable - miniExpectedFrequencies).^2 ./...
            overallExpectedFrequencies([ri,rj],:);
        chi2Values(ri,rj) = sum(residualTable(:));
        chi2Values(rj,ri) = chi2Values(ri,rj);
        
    end
end

%============================


%============================
% GET P-VALUES
%============================

% df = (rows-1)(cols-1), but rows==2 for each comparison
degreesOfFreedom = cols-1;
pValue = 1-chi2cdf(chi2Values,degreesOfFreedom);

