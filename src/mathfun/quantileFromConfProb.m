function q=quantileFromConfProb(p)
% quantileFromConfProb returns the quantile associated to a given confidence probability (for a normal distribution with mean 0 and std 1)
%
% SYNOPSIS      p=quantileFromConfProb(q,m)
%
% INPUT         p               : confidence probability (e.g. 0.95)
%
% OUTPUT        q               : quantile (e.g. 1.96)
%
%
% DEPENDENCES   quantileFromConfProb uses { }
%               quantileFromConfProb is used by {  }
%
% Aaron Ponti, May 11th, 2004

% Check input parameters
if nargin~=1
    error('One input parameter expected');
end

% Value from the inverse cdf for the normal distribution with mean 0 and sigma 1
q=norminv(p+(1-p)/2,0,1);
