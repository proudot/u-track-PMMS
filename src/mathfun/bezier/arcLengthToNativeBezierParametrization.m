function t = arcLengthToNativeBezierParametrization(cP,s,warn)
% function t = arcLengthToNativeBezierParametrization(cP,s)
% arcLengthToNativeBezierParametrization transforms an arc length
% parametrization into the 'native' parametrization of a Bezier curve by
% iteratively splitting the search interval.
%
% Required Inputs:
% cP             A n+1 x M array representing the M-dimensional control
%                points of the Bezier curve of degree n.
%
% s              Any array representing the arc length nodes for which the
%                corresponding 'native' nodes should be computed.          
%
% Optional Inputs:
% warn           If set to 'off' warnings will be suppressed 
%
% Outputs:
% t              Array representing the transformed ('native') nodes
%
% Pascal Berard, December 2011

% Define parameters
iterMax = 1000;
precision = 0.0001;

% Transform line vectors and matrices into column vectors
sizeS = size(s);
s = reshape(s,numel(s),1);

ub = max(1,s); % Upper bound
lb = min(0,s); % Lower bound

if (any(lb < 0) || any(ub > 1)) && (nargin == 2 || strcmp(warn,'on'))
    disp('WARNING: The arc length parameter s is outside the interval [0,1]! Support is experimental!')
end

% Test if all s are between 0 and 1
nS = numel(s); % Number of input nodes
tLow = lb;
tUpper = ub;
sLow = lb;
sUpper = ub;
t = zeros(nS,1); % Output

length = lengthBezier(cP);

for iter=1:iterMax
    tGuess = (tLow+tUpper)/2;
    sGuess = arrayfun(@(a) sign(a)*lengthBezier(cP,0,a)/length,tGuess);
    
    % Determine bounds
    guessTooHigh = sGuess > s;
    tUpper(guessTooHigh) = tGuess(guessTooHigh);
    sUpper(guessTooHigh) = sGuess(guessTooHigh);
    tLow(~guessTooHigh) = tGuess(~guessTooHigh);
    sLow(~guessTooHigh) = sGuess(~guessTooHigh);
    
    % Precision test
    sDelta = arrayfun(@(a,b) lengthBezier(cP,a,b)/length,tLow,tUpper);
    
    ok = sDelta < precision; % Precision reached
    
    % Exit when at least the required precision is reached
    if all(ok)
        % Linear interpolation of t
        f = (s-sLow)/(sUpper-sLow);
        t = tLow + f*(tUpper-tLow);
        break;
    end
end

if iter == iterMax
    disp('WARNING: arcLengthToNativeBezierParametrization: iterMax reached! Result may be inaccurate!');
end

% Restore the input vector/matrix shape
t = reshape(t,sizeS);

end
