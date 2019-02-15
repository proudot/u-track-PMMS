%[h, p, s] = multivarPermTest(X, Y, varargin) performs a multivariate permutation test based on euclidean distance
%
% Inputs:
%              X : n1 x d matrix of samples; d is the dimension.
%              Y : n2 x d matrix of samples
%                  Up to 3-dimensional data are supported in the current implementation.
%
% Options:
%          alpha : Significance level
%           nrep : # of permutations, default: 1900. This gives a coefficient of
%                  variation <=0.10 for alpha = 0.05. Calculated as
%                  nrep = (1-alpha)/cv^2/alpha. See [1] for details.
%  'DistanceFct' : {'euclidan'}|'crossmatch' selects the metric used for the test.
%                  See [2] and [3] for a descritption of each test.
%      'Display' : true|{false} displays the test distribution.
%
% Outputs:
%              h : Selected hypothesis. h = 0 indicates that the null hypothesis cannot be
%                  rejected at significance level 'alpha'.
%              p : p-value
%              s : test statistic
%
% Notes: The crossmatch test is implemented for testing purposes. It is an exact test,
%        the function crossmatchTest() should be used instead.
%
% [1] Efron, B. and Tibshirani, R., "An introduction to the Bootstrap," Ch. 15, 1993.
% [2] Szekely and Rizzo, "Testing for equal distributions in high dimension," InterStat 5, 2004.
% [3] Rosenbaum, "An exact distribution-free test comparing two multivariate
%     distributions based on adjacency," J. R. Statist. Soc. B 67, Part 4, pp. 515-530, 2005.

% Francois Aguet, 08/07/2013

function [hval, pval, s] = multivarPermTest(X, Y, varargin)

ip = inputParser;
ip.addRequired('X', @isnumeric);
ip.addRequired('Y', @isnumeric);
ip.addOptional('alpha', 0.05, @isscalar);
ip.addOptional('nrep', 1900, @isscalar);
ip.addParamValue('DistanceFct', 'euclidean', @(x) any(strcmpi(x, {'euclidean', 'crossmatch'})));
ip.addParamValue('Display', false, @islogical);
ip.parse(X, Y, varargin{:})
nrep = ip.Results.nrep;

n1 = size(X,1);
n2 = size(Y,1);
N = n1+n2;

V = [X; Y];

idx1 = 1:n1;
idx2 = n1+(1:n2);
D = createDistanceMatrix(V, V);

f0 = n1*n2/(n1+n2);

dv = zeros(1,nrep);
switch lower(ip.Results.DistanceFct)
    case 'euclidean'
        % reference distance
        s = f0 * (2/(n1*n2)*sum(sum(D(idx1,idx2)))...
            - 1/n1^2*sum(sum(D(idx1,idx1)))...
            - 1/n2^2*sum(sum(D(idx2,idx2))));
        
        % randomize
        for i = 1:nrep
            idx = randperm(N);
            idx1 = idx(1:n1);
            idx2 = idx(n1+1:N);
            dv(i) = f0 * (2/(n1*n2)*sum(sum(D(idx1,idx2)))...
                - 1/n1^2*sum(sum(D(idx1,idx1)))...
                - 1/n2^2*sum(sum(D(idx2,idx2))));
        end
        pval = sum(dv >= s)/nrep;
    case 'crossmatch' % see crossmatchTest.m
        
        % distance matrix
        [I, J] = meshgrid(1:N);
        E = [nonzeros(tril(I,-1)) nonzeros(tril(J,-1))];
        W = nonzeros(tril(D,-1));
        
        % minimum distance non-bipartite matching using maximum weighted matching
        % -> reverse cost
        W = max(W(:))-W;
        
        % assign all points to unique pairs
        M = maxWeightedMatching(N, E, W);
        % list of pairs:
        m = E(M,:);
        
        % if N is odd, one point remains unassigned -> discard
        if mod(N,2)
            rmIdx = setdiff(1:N, m(:));
            idx = m>rmIdx;
            m(idx) = m(idx)-1;
            if rmIdx<=n1
                n1 = n1-1;
            end
            N = N-1;
        end
        
        % assign labels to pooled samples
        labels = zeros(N,1);
        labels(n1+1:end) = 1;
        
        % test statistic
        t = sum(labels(m),2);
        s = sum(t==1); % # pairs in X,Y
        
        % randomize labels
        for i = 1:nrep
            idx = randperm(N);
            labels = zeros(N,1);
            labels(idx>n1) = 1;
            t = sum(labels(m),2);
            dv(i) = sum(t==1);
        end
        pval = sum(dv <= s)/nrep;
    case 'kernel'
        % not yet implemented
end

hval = pval <= ip.Results.alpha;

if ip.Results.Display
    figure;
    xi = 0:max(dv);
    ni = hist(dv, xi)/numel(dv);
    bar(xi, ni);
    hold on;
    YLim = get(gca, 'YLim');
    plot(s*[1 1], YLim, 'r');
end
