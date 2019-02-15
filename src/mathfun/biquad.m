function f = biquad(x1,x2,a)
%BIQUAD computes the values of a biquadratic function at x1,x2-coordinates
% returns a matrix 
% f(i,j) = a(1) + a(2)*x1(i,j) + a(3)*x2(i,j) + 
%          a(4)*x1(i,j)^2 + a(5)*x2(i,j)^2 +
%          a(6)*x1(i,j)*x2(i,j)
%
% SYNOPSIS s = biquad(x1,x2,a)
%
% INPUT   x1 : matrix with x1 coordinates
%         x2 : matrix with x2 coordinates
%         a  : column vector with coefficients (6 required)
%
% OUTPUT  f  : function values
%          
% SEE ALSO biquadfit

% check the number of coefficients
dima = size(a);
if(dima(1) < dima(2))
   a = permute(a,[2,1]);
end;
[nCoeffs,dummy] = size(a);
if((nCoeffs~=6) || (dummy~=1))
   error('coefficient vector does not contain 6 values');
end;

f = a(1,1)*ones(size(x1)) + ...
   a(2,1)*x1 + a(3,1)*x2 + ...
   a(4,1)*x1.^2 + a(5,1)*x2.^2 + a(6,1)*x1.*x2;

return;