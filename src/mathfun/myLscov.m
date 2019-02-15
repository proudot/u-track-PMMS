function [x,dx,mse] = myLscov(A,b,V)
%MYLSCOV Least squares with known covariance. - returns mse^2 (=sigma0hat^2) as additional output
%   X = MYLSCOV(A,B,V) returns the vector X that solves A*X = B + E
%   where E is normally distributed with zero mean and coVARIANCE V.
%   A must be M-by-N where M > N.  This is the over-determined
%   least squares problem with covariance V of the data. If V is diagonal,
%   the vector containing the diagonal only can be supplied instead of the
%   full matrix.
%   The solution is found without inverting V.
%   
%   Fitting a noise-free constant leads to imaginary values returned for dx
%   and mse. To avoid this, the function returns only the real part of dx
%   and mse. It will return a warning (that can be turned off), if there is
%   a perfect fit.
%
%   The vector X minimizes (A*X-B)'*inv(V)*(A*X-B). The classical
%   linear algebra solution to this problem is:
%
%       X = inv(A'*inv(V)*A)*A'*inv(V)*B
%
%   [X, DX,MSE] =MYLSCOV(A,B,V) returns the standard errors of X in DX. 
%   The standard statistical formula for the standard error of the
%   coefficients is:
%
%      mse = B'*(inv(V) - inv(V)*A*inv(A'*inv(V)*A)*A'*inv(V))*B./(m-n) 
%      dx = sqrt(diag(inv(A'*inv(V)*A)*mse))
%
%   To get the a priori variance (diagonal elements of Q, diag(inv(A'PA)),
%   use dx^2/mse)
%
%   See also SLASH, LSQNONNEG, QR, or LSCOV for Matlab7

%   References:
%       G. Strang, "Introduction to Applied Mathematics",
%       Wellesley-Cambridge, p. 398, 1986.
%
%       F. Graybill, "Theory and Application of the Linear Model",
%       Duxbury Press p. 207, 1976.

%   L. Shure 3-31-89
%   Brad Jones 4-7-95
%   L. Shure 2-8-96
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.15 $  $Date: 2002/04/08 23:51:49 $
%   changed by jonas
%   9/05: made into wrapper for the Matlab7 lscov


% if we're in Matlab 7: Use the built-in lscov
ver=version;
if str2double(ver(1))>=7
    [m,n] = size(A);
    sv = size(V);
    if ~isequal(sv,[m m]),
        if any((sv == m)) && any(sv==1)
            % vector input for lscov7 assumes weights as inverse V
            V = 1./V(:);
        else
            error(sprintf('V must be a %d-by-%d matrix.',m,m));
        end
    end
    [x,dx,mse] = lscov(A,b,V);
    return
end

% Matlab 6

[m,n] = size(A);

% Error checking
if m < n, error('Problem must be over-determined so that M > N.'); end
if m == n, warning('MYLSCOV:NOFIT','There are as many equations as unknowns. There will be no fitting.'),return,end
if size(b,1)~=m, 
  error(sprintf('B must be a column vector of size %d-by-1.',m));
end
sv = size(V);
if ~isequal(sv,[m m]), 
    if any((sv == m))
        V = diag(V);
    else
  error(sprintf('V must be a %d-by-%d matrix.',m,m));
    end
end

[qnull,r] = qr(A);
q = qnull(:,1:n);
r = r(1:n,1:n);
qnull(:,1:n) = [];

g = qnull'*V*qnull;
f = q'*V*qnull;

c = q'*b;
d = qnull'*b;

x = r\(c-f*(g\d));

if nargout > 1
    % Standard formulas are:
    %   mse = b'*(inv(V) - inv(V)*A*inv(A'*inv(V)*A)*A'*inv(V))*b./(m-n); 
    %   dx = sqrt(diag(inv(A'*inv(V)*A)*mse));
    U = (chol(V))';
    z = U\b;
    W = U\A;
    mse = (z'*z - x'*W'*z)/(m-n);    
    [Q, R] = qr(W,0);
    ri = (R\eye(n))';
    dx = (sqrt(sum(ri.*ri)*mse))';
    
    % if all B are equal, there could be imaginary parts of the solution.
    % set real.
    mse = (real(mse));
    dx  = (real(dx));
    
    if mse == 0 | dx == 0
        warning('MYLSCOV:perfectFit','MYLSCOV returned an error-free fit. Subsequent statistical analyses might be wrong')
    end

end
