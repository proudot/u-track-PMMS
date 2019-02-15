%[hval, pval] = ttest2Parameters(mu, s, n, alpha) calculates the two-sample t-test from the input parameters
% Assumes unequal variances, two-tailed test

% Francois Aguet, 2012

function [hval, pval] = ttest2Parameters(mu, s, n, alpha)

if nargin<4;
    alpha = 0.05;
end

v1 = s(1)^2/n(1);
v2 = s(2)^2/n(2);
s_X1X2 = sqrt(v1 + v2);

% test statistic
t = (mu(1)-mu(2)) / s_X1X2;

% degrees of freedom
df = (v1+v2)^2 / (v1^2/(n(1)-1) + v2^2/(n(2)-1));

% two-tailed test
pval = 2 * tcdf(-abs(t), df);

if (pval <= alpha)
    hval = 1;
else
    hval = 0;
end
