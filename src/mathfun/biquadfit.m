function [a , sa ] = biquadfit(x1,x2,f)
%BIQUADFIT least square fit of a biquadratic function
% returns LS estimate of a where a is given by
% f(i,j) + e(i,j) = a(1) + a(2)*x1(i,j) + a(3)*x2(i,j) + 
%                   a(4)*x1(i,j)^2 + a(5)*x2(i,j)^2 +
%                   a(6)*x1(i,j)*x2(i,j)
%
% SYNOPSIS [a , sa ] = biquadfit(x1,x2,f)
%
% INPUT   x1 : matrix with x1 coordinates
%         x2 : matrix with x2 coordinates
%         f  : matrix with function values
%
% OUTPUT  a  : LS estimate of polynom parameters
%         sa : standard errors of the coefficients in a
%          
% SEE ALSO lscov

% check dimensions 
dimX1 = size(x1);
dimX2 = size(x2);
dimF  = size(f);

if((dimX1(1) ~= dimF(1)) || dimX2(1) ~= dimF(1))
   error('dimensions of input matrices not equal');
end;
if((dimX1(2) ~= dimF(2)) || dimX2(2) ~= dimF(2))
   error('dimensions of input matrices not equal');
end;

% store data inb intermediate vectors
vX1 = x1(:);
vX2 = x2(:);
rhs = f(:);

% create design matrix
design = ones(size(rhs));
design = cat(2,design, vX1, vX2, vX1.^2, vX2.^2, vX1.*vX2);

nObs=size(design);
[a,sa] = lscov(design,rhs,eye(nObs));
return;