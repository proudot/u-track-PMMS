function [ rho, theta, X, Y ] = randcirc( n,radius, center )
%randcirc Creates a uniform random distribution of points in a circle of 
% specified radius 
% 
% INPUT (all optional)
% n      : number of points (default: 1)
% radius : radius of the circle (default: 1)
% center : center of the circle (default: [0 0])
%
% OUTPUT
% rho : radius of the points as a column vector
% theta : angle of the points as a column vector
% X : x coordinate of the points as a column vector
% Y : y coordinate of the points as a column vector
%
% EXAMPLE
%
% numPts = 1e5;
% r = 1;
% area = pi*r^2;
% density = numPts / area;
% [~,~,X,Y] = randcirc(numPts,r);
% % plot Ripley's L function
% dist = hypot(X,Y);
% [N,edges] = histcounts(dist,0:r/log(numPts):r)
% plot(edges(2:end),sqrt(cumsum(N)/density/pi) - edges(2:end))
% 

% Mark Kittisopikul
% Jaqaman Lab
% UT Southwestern
% August 21, 2015

if(nargin < 1)
    n = 1;
end
if(nargin < 3)
    center = [ 0 0];
end

% Triangle distribution from 0 to 2 with maximum pdf at 1
% NB: The radially density would otherwise scale with rho if rho were
% simply drawn from a uniform distribution
rho = sum(rand(n,2),2);
% wrap (1,2] to [0,1)
outside = rho > 1;
rho(outside) = 2 - rho(outside);

% Scale by cicle radius
if(nargin > 1)
    % radius is not 1
    rho = rho*radius;
end

% Theta should be a uniform distribution between 0 and 2*pi
if(nargout > 1)
    theta = rand(n,1)*2*pi;
end

% map to cartesian coordinates
if(nargout > 2)
    [X,Y] = pol2cart(theta,rho);
    if(nargin > 2)
        X = X+center(1);
        Y = Y+center(end);
    end
end

end

