%[H T cval] = adtestk(samples, varargin) performs the k-sample Anderson-Darling test
%
% Inputs:
%
%    samples : cell array of sample vectors
%    'Alpha' : alpha value, must be any of [0.25 0.1 0.05 0.025 0.01]. Default: 0.05
%
% Outputs:
%          H : Result of the hypothesis test
%              0: Do not reject H0 at given significance level
%              1: Reject H0 at given significance level
%          T : Value of the test statistic
%       cval : Critical value
%
%
% Implements the test described in
% [1] Scholz and Stephens, J. Am. Stat. Assoc. 82(399), 1987

% Francois Aguet, 03/02/2012

function [H T cval] = adtestkMatlab(samples, varargin)

alphaVec = [0.25 0.1 0.05 0.025 0.01];

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('Alpha', 0.05, @(x) sum(x==alphaVec)==1);
ip.parse(varargin{:});
alpha = ip.Results.Alpha;

k = numel(samples);
n = cellfun(@(i) numel(i), samples)'; % k x 1
Z = sort(horzcat(samples{:}));
Zu = unique(Z);
L = numel(Zu);

N = numel(Z);

H = sum(1./n);
h = sum(1./(1:N-1));


[i,j] = ndgrid(1:N-2,2:N-1);
i = triu(i);
j = triu(j);
g = (N-i).*j;
g(g==0) = [];
g = sum(1./g);

a = (4.0*g-6.0)*(k-1.0) + (10.0-6.0*g)*H;
b = (2.0*g-4.0)*k^2 + 8.0*h*k + (2.0*g-14.0*h-4.0)*H - 8.0*h + 4.0*g - 6.0;
c = (6.0*h+2.0*g-2.0)*k^2 + (4.0*h-4.0*g+6.0)*k + (2.0*h-6.0)*H + 4.0*h;
d = (2.0*h+6.0)*k^2 - 4.0*h*k;


f = zeros(k, L);
for i = 1:k
    [rep, sortedSamples] = getMultiplicity(samples{i});
    f(i, ismember(Zu, sortedSamples)) = rep;
end
% multiplicity of Zu
l = sum(f,1); % 1 x L

% see after Eq. 7
Ma = cumsum(f, 2)-f/2; % k x L
Ba = cumsum(l, 2)-l/2; % 1 x L

% expand all arrays to k x L
Ba = repmat(Ba, [k 1]);
l = repmat(l, [k 1]);
n = repmat(n, [1 L]);

% eq. 7
A2 = (N-1)/N^2 * sum(sum( l./n.*(N*Ma-n.*Ba).^2./(Ba.*(N-Ba)-N*l/4) ));

% variance: depends only on sample sizes (sigma_N^2)
varA = (((a*N + b)*N + c)*N + d) / ( (N-1)*(N-2)*(N-3) );

% test statistic
T = (A2 - (k-1)) / sqrt(varA); % p. 921

% Table I
m = [1 2 3 4 6 8 10];
aIdx = alpha==alphaVec;
if any(m==k-1)
    cTable = [0.326 1.225 1.960 2.719 3.752;
        0.449 1.309 1.945 2.576 3.414;
        0.498 1.324 1.915 2.493 3.246;
        0.525 1.329 1.894 2.438 3.139;
        0.557 1.332 1.859 2.365 3.005;
        0.576 1.330 1.839 2.318 2.920;
        0.590 1.329 1.823 2.284 2.862];
    cval = cTable(m==k-1,aIdx);
else % interpolate critical value
    b = [0.675 1.281 1.645 1.960 2.326;
        -0.245 0.250 0.678 1.149 1.822;
        -0.105 -0.305 -0.362 -0.391 -0.396];
    b = b(:,aIdx);
    cval = b(1) + b(2)/sqrt(k-1) + b(3)/(k-1);
end
H = T >= cval;
