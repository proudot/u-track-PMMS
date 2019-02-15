function [rhoAlpha,slopeIntcpt,dtls]=lineRegression(data,sData,rhoAlpha0,sApprox,sigma)
%LINEREGRESSION computes the model parameters of a line according to total LS
%
% SYNOPSIS [rhoAlpha,slopeIntcpt,dtls]=lineRegression(data,sData,rhoAlpha0,sApprox,sigma)
%
% INPUT data : 2xn matrix with the point coordinates
%       sData: (opt)
%                     [sData_1, sData_2] to define one global uncertainty 
%                     2xn matrix with an uncertainty value for all data points
%       rhoAlpha0: (opt) [rho0, alpha0]
%       sApprox :  (opt) uncertainty of approximations if given;
%                  if rhoAlpha0 ~= [] but sApprox == [], then sApprox is set 
%                  to a large value.
%       sigma :    (opt) stdev of unit weight obervation. Default: 1
%
% OUTPUT  rhoAlpha : estimate of the line parameters
%         slopeIntcpt : parameters tranform3ed to slope and intercept with x_1 = 0;
%         dtls : data structure with the field
%            *.sigma  : aposteriori stdev of unit weight observation
%            *.sRhoAlpha : precision of the parameters propagated through TLS
%            *.sSlopeIntcpt : precision of slope and intercept propagated
%                             from sRhoAlpha
%            *.resid  : 2xn matrix with coordinate residuals
%

% c: Gaudenz Danuser
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=======================
% TEST INPUT
%=======================

if nargin == 0 || isempty(data)
    error('LINEARREGRESSION needs a non-empty input argument ''data''!')
end
[nRows, nPts] = size(data);
% make sure that the data is given as a 2-by-n matrix
if nPts == 2 && nRows ~= 2
    error('input argument ''data'' has to be a 2-by-n array!')
elseif nPts < 3
    error('you need at least 3 points for a linear regression!')
end
    
if nargin < 5 || isempty(sigma)
    sigma = 1;
end

if nargin < 2 || (isempty(sData))
   sData = ones([2,nPts])*sigma;
else
   if(size(sData,1)*size(sData,2) == 2)
       % fill sigma-array
       sData = repmat([sData(1);sData(2)],[1,nPts]);
   end;
end;

if nargin < 3 || (isempty(rhoAlpha0))
   rhoAlpha0 = getApprox(data);
   sApprox = [10^32,10^32];
else
   if nargin < 4 || (isempty(sApprox))
      sApprox = [10^32,10^32];
   end;
end;



[dtls.resid,rhoAlpha,dtls.sRhoAlpha,dtls.sigma,dtls.nIter] = ...
   mexLineRegression(data',sData',rhoAlpha0,sApprox,sigma);

if(rhoAlpha(2) ~= 0)
   sinAlpha = sin(rhoAlpha(2));
   slopeIntcpt(1)=-cos(rhoAlpha(2))/sinAlpha;
   slopeIntcpt(2)=rhoAlpha(1)/sinAlpha;
   % calculate sSlopeIntcpt via error propagation from sRhoAlpha
   sRho = dtls.sRhoAlpha(1);
   sAlpha = dtls.sRhoAlpha(2);
   dtls.sSlopeIntcpt(1) = (sAlpha/sinAlpha^2);
   dtls.sSlopeIntcpt(2) = ...
       sqrt((sRho/sinAlpha)^2 + (rhoAlpha(1)*sAlpha/sinAlpha^2)^2);
else
   slopeIntcpt = [];
end;

dtls.resid = dtls.resid';

%--------------------------------------------------------------------------
function [ra] = getApprox(data)
% determines the parameters over the eigenvalues of the second momentum

m = mean(data,2);
mat = zeros(2);
for(i = 1:size(data,2))
   mat = mat + (data(:,i)-m) * (data(:,i) - m)';
end;

[V,D] = eig(mat);

if(D(1,1) < D(2,2))
   ra(2) = mod(atan2(V(2,1),V(1,1))+pi,pi);
else
   ra(2) = mod(atan2(V(2,2),V(1,2))+pi,pi);
end;
ra(1) = m' * [cos(ra(2));sin(ra(2))];

    
   