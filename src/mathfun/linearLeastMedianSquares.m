function [x,Qxx,goodRows,sigmaB,badRows] = linearLeastMedianSquares(A,B,V,x0,k)
%LINEARLEASTMEDIANSQUARES uses a least median square estimate to do robust least squares fitting of data with outliers
%
% linearLeastMedianSquares fits linear problems in the form of A*x = B + E,
% where E is a hopefully gaussian error term. It first checks for outliers
% in the data (B), and discards those, then it fits the reduced data with a
% standard least squares algorithm (myLscov).
%
% It is possible that removing outliers will turn a column of A into all
% zeros. In this case, the corresponding unknown will be returned as NaN.
%
% see: see Danuser, 1992 or Rousseeuw & Leroy, 1987
%
% SYNOPSIS : [x,Qxx,goodRows,sigmaB,badRows] = linearLeastMedianSquares(A,B,V,x0)
%
% INPUT    : A,B : matrices to describe the least squares problem in the
%                  form A*x = B
%            V   : optional covariance matrix of the values in B (DO NOT INVERT)
%                  You can also supply a vector of covariances, if V is
%                  diagonal.
%                  Note that to calculate leastMedianSquare, only the
%                  diagonal elements of the inverse of V are used!
%            x0  : optional initial guess (default: ones(size(x)))
%                  (fminsearch is sensitive to initial conditions, so it helps
%                  to have a good first guess)
%            k   : outlier threshold. Default is {3} sigma.
%
% OUTPUT   : x,Qxx : fitted unknowns and uncertainty
%            goodRows : rows in B that are not outliers
%            sigmaB   : estimate for the std of the error E without
%                       outliers
%
% c: 3/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%==============
% TEST INPUT
%==============

if nargin < 2 || isempty(A) || isempty(B)
    error('Not enough or empty input arguments in linearLeastMedianSquares');
end

% size A: necessary for defaults
sizA = size(A);

% V: default: eye. However, myLscov accepts a vector, too.
if nargin < 3 || isempty(V)
    V = ones(sizA(1),1);
    diagInvV = ones(sizA(1),1);
    vIsVector = 1;
else
    % if it's a vector, make a diagonal matrix
    sizV = size(V);
    if sizV(1) == sizV(2) && ndims(V) == 2 && sizV(1) == sizA(1)
        % covariance matrix supplied. For median, we only use the diagonal
        diagInvV = diag(inv(V));
        vIsVector = 0;
    elseif any(sizV == sizA(1)) && any(sizV == 1)
        % covariance vector - we can supply that directly into myLscov
        diagInvV = 1./V;
        vIsVector = 1;
    else
        error('incorrect covariance matrix V supplied')
    end
end

% set sigma cutoff
if nargin < 5 || isempty(k)
    k = 3;
end

% --- remove NaNs

% find them
[badRowA,c] = find(~isfinite(A));
[badRowB,c] = find(~isfinite(B));
[badRowV,c] = find(~isfinite(V));
if ~isempty(badRowA) || ~isempty(badRowB)
    badRows = union(badRowA, badRowB);
else
    badRows = [];
end
if ~isempty(badRowV)
    badRows = union(badRows, badRowV);
end

% remove bad rows if there are any
if ~isempty(badRows)

    % make a list of old rows, so that we return the correct goodRows!
    oldRows = [1:sizA(1)];
    oldRows(badRows) = [];

    % remove bad rows
    A(badRows,:) = [];
    B(badRows,:) = [];
    if vIsVector
        V = V(oldRows);
    else
        V = V(oldRows, oldRows);
    end
    diagInvV(badRows) = [];

    % update size of A
    sizA = size(A);
end

%===== END TEST INPUT ======

%===================================
% ROBUST FITTING
%===================================

% If no good initialization provided, we need to fit to find something ok
if nargin < 4 || isempty(x0)
    x0 = fitRoutine(A,B,V,diagInvV,vIsVector,ones(sizA(2),1),k);
else
    if length(x0)~=sizA(2)
        error('not the right number of initial guesses!')
    end
end

% fit. If there were NaNs in the initial guess, set them to 0 - otherwise
% we have a problem.
x0(isnan(x0)) = 0;
[x,Qxx,goodRows,sigmaB,badR] = fitRoutine(A,B,V,diagInvV,vIsVector,x0,k);

%====================
% OUTPUT
%====================

%set the correct goodRows
if ~isempty(badRows)
    goodRows = oldRows(goodRows);
    badR = oldRows(badR);
    badRows = union(badR,badRows);
else
    badRows = badR;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== FITTING FUNCTION =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,Qxx,goodRows,sigmaB,badRows] = fitRoutine(A,B,V,diagInvV,vIsVector,x0,k)

%========================
% LEAST MEDIAN SQUARES
%========================

% define magic number:
magicNumber2=1.4826^2; %value important for calculation of sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987

% generate function string
functionString = ['inline(''median((A*x-B).*diagInvV.*(A*x-B))'',''x'',''A'',''B'',''diagInvV'')'];
% generate rest of fminsearch call
fminCall = ['fminsearch(' functionString ',x0,options,A,B,diagInvV)'];

% set options
options = optimset('Display','off');

% minimize
xFmin = eval(fminCall);

% calculate statistics
res2 = (A*xFmin-B).^2;
medRes2 = nanmedian(res2);

%testvalue to calculate weights
testValue=res2/(magicNumber2*medRes2);

%goodRows: weight 1, badRows: weight 0
goodRows=find(testValue<=k^2);
badRows=find(testValue>k^2);

% ssq=sum(res2);
nGoodRows = length(goodRows);
if nGoodRows > 4
    sigmaB=sqrt(sum(res2(goodRows))/(nGoodRows-4));
else
    sigmaB = NaN;
end

%====END LMS=========


%=======================
% LINEAR LEAST SQUARES
%=======================

if ~isempty(goodRows)

% check whether we have a zero-column in A
zeroCols = all(A(goodRows,:)==0,1);
% remove the zero-columns
A(:,zeroCols) = [];

% call myLscov
if vIsVector
    [x,Qxx,mse] = myLscov(A(goodRows,:),B(goodRows,:),V(goodRows));
else
    [x,Qxx] = myLscov(A(goodRows,:),B(goodRows,:),V(goodRows,goodRows));
end

% if we removed zero-cols, we have to insert NaNs now.
if any(zeroCols)
    xTmp = repmat(NaN,size(xFmin));
    xTmp(~zeroCols) = x;
    x = xTmp;
    % also assign NaN to Qxx where applicable
    xTmp(~zeroCols) = Qxx;
    Qxx = xTmp;
end
else
    % all data is considered outlier (possible if highly bimodal data)
    [x,Qxx] = deal([]);
end

%=== END LSQ ========

