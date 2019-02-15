%[h, p, s] = crossmatchTest(X, Y, varargin) implements the cross-match test for comparing two multi-variate distributions
%
% Inputs:
%         X : n1 x d matrix of samples; d is the dimension.
%         Y : n2 x d matrix of samples
%             Up to 3-dimensional data are supported in the current implementation.
%
% Options:
%     alpha : Significance level
%  'Display': true|{false} for 2D data, display pairwise assignments of samples
%
% Outputs:
%         h : Selected hypothesis. h = 0 indicates that the null hypothesis cannot be
%             rejected at significance level 'alpha'.
%         p : p-value
%         s : value of the test statistic (= number of matched pairs)
%
%
% Implements the cross-match test described in
% P.R. Rosenbaum, "An exact distribution-free test comparing two multivariate
% distributions based on adjacency," J. R. Statist. Soc. B 67, Part 4, pp. 515-530, 2005
%
% This implementation was partially adapted from an implementation for R,
% available at http://cran.r-project.org/web/packages/crossmatch/index.html
%
% Examples:
% 
% [hval, pval] = crossmatchTest(rand(10,2), rand(15,2), 0.05, 'Display', true)
% [hval, pval] = crossmatchTest(rand(10,2), 1+rand(15,2), 0.05, 'Display', true)

% Francois Aguet, 08/10/2013

function [h, p, a1] = crossmatchTest(X, Y, varargin)

ip = inputParser;
ip.addRequired('X');
ip.addRequired('Y');
ip.addOptional('alpha', 0.05, @isscalar);
ip.addParamValue('Display', false, @islogical);
ip.parse(X, Y, varargin{:})

[n1,dim] = size(X);

if size(Y,2)~=dim || dim>4
    error('Inputs X and Y must have the same dimension. Up to 3D data are supported.');
end

n2 = size(Y,1);
N = n1+n2;

% pool observations
v = [X; Y];

% assign labels to pooled input
labels = zeros(N,1);
labels(n1+1:end) = 1;

oddN = mod(N,2);
if mod(N,2) % add point for matching
    v = [v; NaN(1,dim)];
    N = N+1;
end

% distance matrix
D = createDistanceMatrix(v,v);
[I, J] = meshgrid(1:N);
E = [nonzeros(tril(I,-1)) nonzeros(tril(J,-1))];
W = D(tril(ones(N),-1)~=0);

% minimum distance non-bipartite matching using maximum weighted matching
% -> reverse cost
W = max(W(:))-W;
W(isnan(W)) = 0; % no cost to match to placeholder

% assign all points to unique pairs
M = maxWeightedPerfectMatching(N, E, W);
% list of pairs:
m = E(M,:);

% discard dummy pair if odd #points
if oddN
    rmIdx = any(m==N,2);
    rmEl = setdiff(m(rmIdx,:), N);
    m(rmIdx,:) = [];
    N = N-2;
    if rmEl<=n1
        n1 = n1-1;
    else
        n2 = n2-1;
    end
end

% # of pairs
t = sum(labels(m),2);
a0 = sum(t==0); % # pairs in X
a1 = sum(t==1); % # pairs in X,Y
a2 = sum(t==2); % # pairs in Y
I = a0+a1+a2;

% vector of possible X,Y pair counts up to a1
if mod(n1,2)
    av = 1:2:a1;
else
    av = 0:2:a1;
end

if N<340
    tmp = factorial(I)/nchoosek(N,n1);
    A0 = I-(n1+av)/2;
    A2 = (n1-av)/2;
    p = tmp *  2.^av./(factorial(A0).*factorial(av).*factorial(A2));
    p = cumsum(p);
else
    % Distribution under H0: normal distribution with (mu, sigma)
    % Expectation and SD for a1:
    mu = n1*n2/(N-1);
    sigma = sqrt(2*n1*(n1-1)*n2*(n2-1)/((N-3)*(N-1)^2));
    
    % In the paper/R code, p = normcdf(a1, mu, sigma) is used.
    % However, the exact distribution is defined only at even or odd integers,
    % depending on parity of input. The following is the accurate calculation:
    p = 2*cumsum(normpdf(av, mu, sigma));
end
p = p(end);

h = p < ip.Results.alpha;

if ip.Results.Display && dim==2
    figure; hold on;
    plot([v(m(:,1),1) v(m(:,2),1)]', [v(m(:,1),2) v(m(:,2),2)]', 'k');
    plot(X(:,1), X(:,2), 'b.');
    plot(Y(:,1), Y(:,2), 'r.');
end
