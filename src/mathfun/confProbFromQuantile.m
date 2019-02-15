function p=confProbFromQuantile(q)
% confProbFromQuantile returns the confidence probability associated to a given quantile (for a normal distribution with mean 0 and std 1)
%
% SYNOPSIS      p=confProbFromQuantile(q,m)
%
% INPUT         q               : quantile (e.g. 1.96)
%
% OUTPUT        p               : confidence probability (e.g. 0.95)
%
% DEPENDENCES   confProbFromQuantile uses { }
%               confProbFromQuantile is used by {  }
%
% Aaron Ponti, May 11th, 2004

% Check input parameters
if nargin~=1
    error('One input parameter expected');
end

% Densely sampled probability range
probs=[0:0.001:1];

% Inverse cdf for the normal distribution with mean 0 and sigma 1
icdf=norminv(probs,0,1);

% Find the value of the icdf closest to the input quantile
d=abs(icdf-q);
indx=find(d==min(d));

% Return the corresponding confidence probability
p=probs(indx);
p=1-2*(1-p);
