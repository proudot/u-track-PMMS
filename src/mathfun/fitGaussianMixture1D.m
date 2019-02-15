%[mu, sigma, A] = fitGaussianMixture1D(data, n, varargin)
%
% Inputs:
%            data : samples of a distribution
%               n : number of Gaussians to fit
%     {'Display'} : 'on' | {'off'} to view results of the fit
%
% Options:
%     'Optimizer' : 'fmincon'|{'lsqnonlin'}
%                   With fmincon, the amplitudes are constrained to sum(A)==1
%  ConstrainMeans : true|{false} constrains the means to multiples of the first mixture
%                   and the amplitudes to monotonically decrease, i.e., A1>A2, A2>A3,...
%
% Outputs:
%              mu : means of the Gaussians. 2nd row contains propagated error
%           sigma : standard deviations of the Gaussians. 2nd row: propagated error
%               A : amplitudes/relative contributions of the Gaussians. 2nd row: propagated error
%
% Note: when constraining the means, the minimum bound for the lowest mean is 
%       set to the 5th percentile of the data. This value can be adjusted with
%       the parameter 'MinPercentile'. Default 5.

% Francois Aguet, 07/19/2011 (last updated: 07/12/2013)

function [mu, sigma, A, RSS, BIC, prmStd] = fitGaussianMixture1D(arg1, arg2, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('arg1', @isvector);
ip.addRequired('arg2');
ip.addOptional('arg3', []);
ip.addParamValue('Init', [], @isvector);
ip.addParamValue('Display', 'on', @(x) any(strcmpi(x, {'on','off'})));
ip.addParamValue('Optimizer', 'lsqnonlin', @(x) any(strcmpi(x, {'fmincon','lsqnonlin'})));
ip.addParamValue('ConstrainMeans', false, @islogical);
ip.addParamValue('ConstrainSD', false, @islogical);
ip.addParamValue('MinPercentile', 5);
ip.parse(arg1, arg2, varargin{:});
init = ip.Results.Init;

% determine whether inputs were samples (CDF fit) or a density (PDF fit)
if ~isempty(ip.Results.arg3)
    x = ip.Results.arg1;
    f = ip.Results.arg2;
    n = ip.Results.arg3;
    % make sure data is correctly normalized
    dx = x(2)-x(1);
    f = f/sum(f)/dx;
    mode = 'PDF';
    mu0 = sum(x.*f)/sum(f);
    sigma0 = sqrt(sum((x-mu0).^2.*f)/sum(f));
    [cdf,idx] = unique(cumsum(f)*dx);
    Amin = interp1(cdf, x(idx), ip.Results.MinPercentile/100);
else
    data = ip.Results.arg1;
    n = ip.Results.arg2;
    [f, x] = ecdf(data);
    mode = 'CDF';
    mu0 = mean(data);
    sigma0 = std(data);
    Amin = prctile(data, ip.Results.MinPercentile);
end

opts = optimset('MaxFunEvals', 1e4, ...
    'MaxIter', 1e4, ...
    'Display', 'off', ...
    'TolX', 1e-10, ...
    'Tolfun', 1e-10);


if isempty(init)
    mu_init = cumsum(ones(1,n)/(n+1))*2*mu0;
    sigma_init = sigma0/sqrt(n)*ones(1,n);   
    A_init = (n:-1:1)/n/(n+1)*2;
    if ip.Results.ConstrainMeans
        if ~ip.Results.ConstrainSD
            % [mu1, A1, mu2, A2, ..., sigma]
            init = [mu_init(1) reshape([sigma_init; A_init], [1 2*n])];
            % constrain minimum mean to 5th percentile of data
            lb = [Amin zeros(1,2*n)];
            ub = [Inf repmat([Inf 1], [1 n])];
        else
            % mu, sigma, A1, ... An
            init = [mu_init(1) sigma_init(1) A_init];
            lb = [-Inf 0 zeros(1,n)];
            ub = [Inf Inf ones(1,n)]; 
        end
    else
        lb = repmat([-Inf 0 0], [1 n]);
        ub = repmat([Inf Inf 1], [1 n]);
        init = reshape([mu_init; sigma_init; A_init], [1 3*n]);
    end
end

switch mode
    case 'CDF'
        cost = @costCDF;
        costConstr = @costCDFconstr;
        costConstr2 = @costCDFconstr2;
    case 'PDF'
        cost = @costPDF;
        costConstr = @costPDFconstr;
        costConstr2 = @costPDFconstr2;
end

switch ip.Results.Optimizer
    case 'fmincon'
        if ip.Results.ConstrainMeans
            if ~ip.Results.ConstrainSD
                Aeq = [0 repmat([0 1], [1 n])]; % equality constraint: sum(A)==1
                beq = 1;
                if n>1
                    % Components must have diminishing weights: A_{n+1} <= A_{n}
                    Aneq = arrayfun(@(i) circshift([0 -1 0 1 zeros(1,2*(n-2))], [0 2*i]), 0:n-2, 'unif', 0);
                    Aneq = [zeros(n-1,1) vertcat(Aneq{:})];
                    bneq = zeros(n-1,1);
                else
                    Aneq = [];
                    bneq = [];
                end
                [p,RSS] = fmincon(@(i) sum(costConstr(i, x, f).^2), init, Aneq, bneq, Aeq, beq,...
                    lb, ub, [], optimset(opts, 'Algorithm', 'interior-point'));
                mu = p(1)*(1:n);
                sigma = p(2:2:end);
                A = p(3:2:end);
            else % mean and SD constrained
                Aeq = [0 0 ones(1,n)];
                beq = 1;
                if n>1
                    Aneq = arrayfun(@(i) circshift([-1 1 zeros(1,n-2)], [0 i]), 0:n-2, 'unif', 0);
                    Aneq = [zeros(n-1,2) vertcat(Aneq{:})];
                    bneq = zeros(n-1,1);
                else
                    Aneq = [];
                    bneq = [];
                end
                [p,RSS] = fmincon(@(i) sum(costConstr2(i, x, f).^2), init, Aneq, bneq, Aeq, beq,...
                    lb, ub, [], optimset(opts, 'Algorithm', 'interior-point'));
                mu = p(1)*(1:n);
                sigma = p(2)*sqrt(1:n);
                A = p(3:end);
            end
        else
            % Aeq * x = beq
            % equality constraint: sum(A)==1
            Aeq = repmat([0 0 1], [1 n]); 
            beq = 1;
            Aneq = [];
            bneq = [];
            [p,RSS] = fmincon(@(i) sum(cost(i, x, f).^2), init, Aneq, bneq, Aeq, beq, lb, ub, [], optimset(opts, 'Algorithm', 'interior-point'));
            mu = p(1:3:end);
            sigma = p(2:3:end);
            A = p(3:3:end);
        end
    case 'lsqnonlin'
        if ip.Results.ConstrainMeans
            if ~ip.Results.ConstrainSD
                [p,RSS,~,~,~,~,J] = lsqnonlin(costConstr, init, lb, ub, opts, x, f);
                mu = p(1)*(1:n);
                sigma = p(2:2:end);
                A = p(3:2:end);
            else
                [p,RSS,~,~,~,~,J] = lsqnonlin(costConstr2, init, lb, ub, opts, x, f);
                mu = p(1)*(1:n);
                sigma = p(2)*sqrt(1:n);
                A = p(3:end);
            end
        else
            [p,RSS,~,~,~,~,J] = lsqnonlin(cost, init, lb, ub, opts, x, f);
            mu = p(1:3:end);
            sigma = p(2:3:end);
            A = p(3:3:end);
        end
end
ns = numel(f);
BIC = ns*log(RSS/ns) + numel(p)*log(ns);

% covariance matrix, error propagation
C = RSS*full(inv(J'*J));
p_std = sqrt(diag(C)/(ns-length(p) - 1))';
if ip.Results.ConstrainMeans
    prmStd.mu = p_std(1);
    prmStd.sigma = p_std(2:2:end);
    prmStd.A = p_std(3:2:end);
else
    prmStd.mu = p_std(1:3:end);
    prmStd.sigma = p_std(2:3:end);
    prmStd.A = p_std(3:3:end);
end


if strcmpi(ip.Results.Display, 'on')
    if strcmpi(mode, 'CDF')
        figure;
        hold on;
        plot(x, f, 'k', 'LineWidth', 1);
        for i = 1:n
            f = mixtureModelCDF(x, mu(i), sigma(i), A(i));
            plot(x, f, 'b--', 'LineWidth', 1);
        end
        plot(x, mixtureModelCDF(x, mu, sigma, A), 'r--', 'LineWidth', 1.5);
        axis([min(x) max(x) 0 1]);
    end
    
    % always plot PDF
    figure;
    hold on;
    if strcmpi(mode, 'PDF')
        plot(x, f, 'k', 'LineWidth', 1);
    else
        [di,x] = ksdensity(data, 'npoints', 100);
        plot(x, di, 'k', 'LineWidth', 1);
        title('Kernel density');
    end
    for i = 1:n
        plot(x, mixtureModelPDF(x, mu(i), sigma(i), A(i)), 'b--', 'LineWidth', 1.5);
    end
    plot(x, mixtureModelPDF(x, mu, sigma, A), 'r--', 'LineWidth', 1.5);
end


function v = costCDF(p, x, f)
mu = p(1:3:end);
sigma = p(2:3:end);
A = p(3:3:end);
v = mixtureModelCDF(x, mu, sigma, A) - f;

function v = costCDFconstr(p, x, f)
n = (numel(p)-1)/2;
mu = p(1)*(1:n);
sigma = p(2:2:end);
A = p(3:2:end);
v = mixtureModelCDF(x, mu, sigma, A) - f;

function v = costCDFconstr2(p, x, f)
n = numel(p)-2;
mu = p(1)*(1:n);
sigma = p(2)*sqrt(1:n);
A = p(3:end);
v = mixtureModelCDF(x, mu, sigma, A) - f;

function f = mixtureModelCDF(x, mu, sigma, A)
f = zeros(size(x));
for i = 1:numel(A)
    f = f + A(i)*normcdf(x, mu(i), sigma(i));
end


function v = costPDF(p, x, f)
mu = p(1:3:end);
sigma = p(2:3:end);
A = p(3:3:end);
v = mixtureModelPDF(x, mu, sigma, A) - f;

function v = costPDFconstr(p, x, f)
n = (numel(p)-1)/2;
mu = p(1)*(1:n);
sigma = p(2:2:end);
A = p(3:2:end);
A = sort(A, 'descend'); % redundant for fmincon case
v = mixtureModelPDF(x, mu, sigma, A) - f;

function v = costPDFconstr2(p, x, f)
n = numel(p)-2;
mu = p(1)*(1:n);
sigma = p(2)*sqrt(1:n);
A = p(3:end);
A = sort(A, 'descend'); % redundant for fmincon case
v = mixtureModelPDF(x, mu, sigma, A) - f;

function f = mixtureModelPDF(x, mu, sigma, A)
f = zeros(size(x));
for i = 1:numel(A)
    f = f + A(i)*normpdf(x, mu(i), sigma(i));
end
