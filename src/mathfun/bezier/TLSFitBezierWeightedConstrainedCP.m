function [P t res] = TLSFitBezierWeightedConstrainedCP(data, w, n, maxCurvature, varargin)
% TLSFitBezierWeightedConstrainedCP computes the total weighted least 
% squares fit of a nth-degree Bezier curve to a set of data points. 
%
% Required Inputs:
% data           A m x d array representing a set of d-dimensional points
% w              A m x d array representing a set of d-dimensional weights
% n              Degree of the Bezier curve
% maxCurvature   An approximate upper bound on the curvature of the Bezier curve
%
% Optional Inputs:
% MaxFunEvals    Maximum number of fonctional evaluations during lsqnonlin or lsqlin.
% MaxIter        Maximum number of interations during lsqnonlin or lsqlin.
% Display        Verbose mode during lsqnonlin or lsqlin.
% TolX           Tolerance on the solution (i.e. t) during lsqnonlin or lsqlin.
% TolFun         Tolerance on the functional during lsqnonlin or lsqlin.
%
% Outputs:
% P              A n+1 x d array representing the set of d-dimensional
%                control points defining the Bezier curve.
%
% t              A m x 1 array of parameter value for each input points. t
%                belongs to [0, 1].
%
% res            A m x 1 array of residue, i.e the components of the distance
%                between a data point and the fitted curve.
%
% Pascal Berard, August 2011

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isnumeric);
ip.addRequired('n', @(n) n >= 0);
ip.addParamValue('MaxFunEvals', 1e4, @isscalar);
ip.addParamValue('MaxIter', 200, @isscalar);
ip.addParamValue('Display', 'off', @isstring);
ip.addParamValue('TolX', 1e-8, @isscalar);
ip.addParamValue('TolFun', 1e-8, @isscalar);

ip.parse(data, n, varargin{:});
maxFunEvals = ip.Results.MaxFunEvals;
maxIter = ip.Results.MaxIter;
display = ip.Results.Display;
tolX = ip.Results.TolX;
tolFun = ip.Results.TolFun;

% Define options of lsqnonlin algorithm
opts = optimset('Jacobian', 'on', ...
    'MaxFunEvals', maxFunEvals, ...
    'MaxIter', maxIter, ...
    'Display', display, ...
    'TolX', tolX, ...
    'DerivativeCheck', 'off', ...
    'LargeScale', 'off', ...
    'Algorithm', 'levenberg-marquardt', ...
    'Tolfun', tolFun);

[m dim] = size(data);

if n == 0
    % Isotropic Gaussian fit
    d = data.*w;
    [Px, ~, resx] = lsqlin(w(:,1), d(:,1), [], [], [], [], [], [], data(1,1), opts);
    [Py, ~, resy] = lsqlin(w(:,2), d(:,2), [], [], [], [], [], [], data(1,2), opts);
    [Pz, ~, resz] = lsqlin(w(:,3), d(:,3), [], [], [], [], [], [], data(1,3), opts);
    P = [Px,Py,Pz];
    res = [resx,resy,resz]./w;
    t = zeros(m,1);
else
    Cnk = diag([1 cumprod(n:-1:1) ./ cumprod(1:n)]);
    
    % Fit a line
    % TODO: Replace with a weighted TLS fit
    meanW = mean(w,1);
    [x0,a,~,~] = ls3dline(bsxfun(@times,data,meanW));
    [minPoint,maxPoint,~,~,~] = projectPointsOntoLine(data,x0'./meanW,a'./meanW);
    
    % Compute the planar constraints
    planeNormal = (maxPoint-minPoint);
    planePoints = repmat(minPoint,n+1,1) + repmat(1/n*(0:n)',1,dim) .* repmat(planeNormal,n+1,1);
    d = planePoints * planeNormal';
    
    Aeq_row = [planeNormal; zeros(n,dim)];
    Aeq_row = Aeq_row(:)';
    Aeq = zeros(n+1,(n+1)*dim);
    Aeq(1,:) = Aeq_row;
    for k=1:n
        Aeq(k+1,:) = circshift(Aeq_row,[0,k]);
    end
    beq = d;
    
    % Minimal curvature radius
    minRad = 1/maxCurvature;
    
    % Deviation constraints
    len = norm(planeNormal);
    if len < 2*minRad
        alpha = minRad - sqrt(minRad^2-0.25*len^2);
    else
        alpha =  len/2;
    end
    
    ub = planePoints+alpha;
    ub = ub(:);
    lb = planePoints-alpha;
    lb = lb(:);
    
    % TODO: Generalize to n dimensions
    if dim == 2
        [I J] = meshgrid([1,-1],[1,-1]);
        C = [I(:), J(:)];
    elseif dim == 3
        [I J K] = meshgrid([1,-1],[1,-1],[1,-1]);
        C = [I(:), J(:), K(:)];
    end
    
    block = zeros(2^dim,(n+1)*dim);
    block(:,1:n+1:dim*(n+1)) = C;
    A = zeros((n+1)*2^dim,(n+1)*dim);
    A(1:2^dim,:) = block;
    for d=1:n
        A(d*2^dim+1:(d+1)*2^dim,:) = circshift(block,[0,d]);
    end
    
    b = sum(reshape(repmat(planePoints',2^dim,1),dim,2^dim*(n+1))'.*repmat(C,n+1,1),2)+sqrt(2)*alpha;
    
    % Define initial node vector t
    t = linspace(0,1,m)';
    
    % Initial control points P
    P = planePoints;
    
    % Compute the weights diagonal matrix
    W = sparse(diag(w(:)));
    
    % Value of the objective function
    resnormOld = -1;
    
    for i=1:maxIter
        % Solve the non-linear optimization on t
        Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
        fun = @(t) r(t, data, W, Cnk, Cn_1k, n, P);
        lbNonLin = []; ubNonLin = [];
        [t, ~, res] = lsqnonlin(fun, t, lbNonLin, ubNonLin, opts);
        
        % Compute Bernstein Matrix
        B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;
        z = zeros(size(B));
        
        % Solve the linear optimization on P
        if n ~= 1

            % TODO: Generalize to n dimensions
            if dim == 2
                C = W*[B z; z B];
            elseif dim == 3
                C = W*[B z z; z B z; z z B];
            end
            
            d = W*data(:);
            [P, resnorm, res] = lsqlin(C, d, A, b, Aeq, beq, lb, ub, P(:), opts);
        else
            % TODO: Generalize to n dimensions
            if dim == 2
                C = W*[B z; z B];
            elseif dim == 3
                C = W*[B z z; z B z; z z B];
            end
            
            d = W*data(:);
            [P, resnorm, res] = lsqlin(C, d, [], [], [], [], [], [], P(:), opts);
        end
        
        P = reshape(P(1:end), [n+1,dim]);

        if resnorm-resnormOld < tolFun
            break;
        else
            resnormOld = resnorm;
        end
    end
    
    % Solve the non-linear optimization on t
    Cn_1k = diag([1 cumprod(n-1:-1:1) ./ cumprod(1:n-1)]);
    fun = @(t) r(t, data, W, Cnk, Cn_1k, n, P);
    lbNonLin = []; ubNonLin = [];
    [t, ~, res] = lsqnonlin(fun, t, lbNonLin, ubNonLin, opts);
    
    % Truncate the Bezier curve
    [P,t] = truncateBezier(P,min(t),max(t),t);
    
    % Compute unweighted residuals
    res = res./w(:);
    res = reshape(res,[m, dim]);
        
end

function [F J] = r(t, data, W, Cnk, Cn_1k, n, P)

[m dim] = size(data);

% Compute Bernstein Matrix
B = (bsxfun(@power, t, 0:n) .* bsxfun(@power, 1 - t, n:-1:0)) * Cnk;

% Compute the residuals
F = data - B * P;
F = W * F(:);

if nargout > 1
    % Compute the Benstein Matrix of order n-1
    Bn_1 = (bsxfun(@power, t, 0:n-1) .* bsxfun(@power, 1 - t, n-1:-1:0)) * Cn_1k;
    
    % Compute derivative of B against t
    z = zeros(m,1);
    Bt = n * ([z Bn_1] - [Bn_1 z]);
    
    % Compute the Jacobian matrix
    % TODO: Generalize to n dimensions
    if dim == 2
        J11 = W(1:m,1:m)*diag(-Bt * P(:,1));
        J21 = W(m+1:end,m+1:end)*diag(-Bt * P(:,2));
        J = [J11; J21];
    elseif dim == 3
        J11 = W(1:m,1:m)*diag(-Bt * P(:,1));
        J21 = W(m+1:2*m,m+1:2*m)*diag(-Bt * P(:,2));
        J31 = W(2*m+1:end,2*m+1:end)*diag(-Bt * P(:,3));
        J = [J11; J21; J31];
    end
end


