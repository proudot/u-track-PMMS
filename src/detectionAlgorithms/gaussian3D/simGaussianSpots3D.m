% SIMGAUSSIANSPOTS3D generates 3D Gaussian spots in a volume.
%
%   Input:
%           dims:  data dimensions: x, y, z
%           sigma: standard deviation (in pixels). If entered as two-element
%                  vector, sigma(1) is the x-y s.d., and sigma(2) the z s.d.
%
%   Options ('specifier', value pairs):
%
%           'X': Nx3 matrix with center coordinates of 3D Gaussians (can be subpixel)
%           'A': vector of amplitudes, or single scalar value (used for all Gaussians).
%           'npoints'   : number of 2D Gaussian to be generated
%           'Background': value of background offset for entire volume
%           'Border' : border conditions: 'padded' (default), 'periodic', or 'truncated'
%           'Normalization: 'on' | {'off'} divides Gaussians by 2*pi*sigma^2 when 'on'
%           'Verbose': 'on' | {'off'}
%
%   Output:
%           frame: volume with 3D Gaussian signals
%           X:     center coordinates of the Gaussians: [x, y, z]
%           sv:    vector of standard deviations
%           Av:    vector of amplitudes
%
% Example:
% vol = simGaussianSpots3D([200 100 50], 2, 'npoints', 50, 'Border', 'periodic');

% Francois Aguet, 07/2013

function [frame, X, sv, Av] = simGaussianSpots3D(dims, sigma, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('dims', @isnumeric);
ip.addRequired('sigma', @(x) isnumeric(x));
ip.addParamValue('X', [], @(x) isnumeric(x) && size(x,2)==3);
ip.addParamValue('A', []);
ip.addParamValue('npoints', 1);
ip.addParamValue('Background', 0);
ip.addParamValue('Window', []);
ip.addParamValue('Verbose', 'off', @(x) any(strcmpi(x, {'on', 'off'})));
ip.addParamValue('Border', 'padded', @(x) any(strcmpi(x, {'padded', 'periodic', 'truncated'})));
ip.addParamValue('Normalization', 'off', @(x) any(strcmpi(x, {'analytical', 'sum', 'off'})));
ip.addParamValue('NonOverlapping', false, @islogical);
ip.parse(dims, sigma, varargin{:});

np = ip.Results.npoints;
c = ip.Results.Background;

if ~isempty(ip.Results.X)
    xv = ip.Results.X(:,1);
    yv = ip.Results.X(:,2);
    zv = ip.Results.X(:,3);
else
    xv = [];
    yv = [];
    zv = [];
end

if ~isempty(xv)
    np = size(xv,1);
end

Av = ip.Results.A(:);
S = ip.Results.sigma(:);

% feature vectors
if numel(S)==1
    sv = S*ones(np,1);
    rv = S*ones(np,1);
elseif numel(S)==2
    sv = S(1)*ones(np,1);
    rv = S(2)*ones(np,1);
else
    sv = S(:,1);
    rv = S(:,2);
end
nx = dims(1);
ny = dims(2);
nz = dims(3);

% window size for each source
wv = ip.Results.Window;
if isempty(wv)
    wv = ceil([4*sv 4*rv]);
    % elseif numel(wv)==1
    %     wv = wv*ones(size(sv));
end
sf = wv(1,1) ./ wv(1,2); % scaling between x,z assumed identical for all points

if isempty(Av)
    Av = ones(np,1);
elseif numel(Av)==1
    Av = Av*ones(np,1);
end

% Generate point coordinates, if input is empty
if strcmpi(ip.Results.Border, 'padded');
    % discard signals close to image border
    if ~isempty(xv)
        idx = xv <= wv(:,1) | yv <= wv(:,1) | zv <= wv(:,2) |...
            xv > nx-wv(:,1) | yv > ny-wv(:,1) | zv > nz-wv(:,2);
        xv(idx) = [];
        yv(idx) = [];
        zv(idx) = [];
        Av(idx) = [];
        sv(idx) = [];
        rv(idx) = [];
        if strcmpi(ip.Results.Verbose, 'on')
            fprintf('Number of discarded points: %d\n', numel(idx));
        end
        np = length(xv);
    else
        if ip.Results.NonOverlapping
            xv = [];
            yv = [];
            zv = [];
            while numel(xv)<np
                xcand = (nx-2*wv(:,1)-1).*rand(np,1) + wv(:,1)+1;
                ycand = (ny-2*wv(:,1)-1).*rand(np,1) + wv(:,1)+1;
                zcand = (nz-2*wv(:,2)-1).*rand(np,1) + wv(:,2)+1;

                idx = KDTreeBallQuery([xcand ycand zcand*sf; xv yv zv*sf], [xcand ycand zcand*sf], 1.5*wv(:,1));
                
                idx(cellfun(@numel, idx)>1) = [];
                idx = vertcat(idx{:});
                xv = [xv; xcand(idx)]; %#ok<AGROW>
                yv = [yv; ycand(idx)]; %#ok<AGROW>
                zv = [zv; zcand(idx)]; %#ok<AGROW>
            end
            xv = xv(1:np);
            yv = yv(1:np);
            zv = zv(1:np);
        else
            xv = (nx-2*wv(:,1)-1).*rand(np,1) + wv(:,1)+1;
            yv = (ny-2*wv(:,1)-1).*rand(np,1) + wv(:,1)+1;
            zv = (nz-2*wv(:,2)-1).*rand(np,1) + wv(:,2)+1;
        end
        [~, idx] = sort(xv+yv*nx); % sort spots according to row position
        xv = xv(idx);
        yv = yv(idx);
        zv = zv(idx);
    end
else
    if ~isempty(xv) && (any(min([xv; yv; zv])<0.5) || any(max(xv)>=nx+0.5) || any(max(yv)>=ny+0.5) || any(max(zv)>=nz+0.5))
        error('All points must lie within x:[0.5 nx+0.5), y:[0.5 ny+0.5), z:[0.5 nz+0.5].');
    end
    if isempty(xv)
        if ip.Results.NonOverlapping
            xv = [];
            yv = [];
            zv = [];
            while numel(xv)<np
                xcand = nx*rand(np,1)+0.5;
                ycand = ny*rand(np,1)+0.5;
                zcand = nz*rand(np,1)+0.5;
                idx = KDTreeBallQuery([xcand ycand zcand; xv yv zv], [xcand ycand zcand], 1.5*max(wv,[],2));
                idx(cellfun(@numel, idx)>1) = [];
                idx = vertcat(idx{:});
                xv = [xv; xcand(idx)]; %#ok<AGROW>
                yv = [yv; ycand(idx)]; %#ok<AGROW>
                zv = [zv; zcand(idx)]; %#ok<AGROW>
            end
            xv = xv(1:np);
            yv = yv(1:np);
            zv = zv(1:np);
        else
            xv = nx*rand(np,1)+0.5;
            yv = ny*rand(np,1)+0.5;
            zv = nz*rand(np,1)+0.5;
        end
    end
end

% background image
frame = c*ones(ny, nx, nz);

xi = round(xv);
yi = round(yv);
zi = round(zv);
dx = xv-xi;
dy = yv-yi;
dz = zv-zi;

X = [xv(:) yv(:) zv(:)];

if strcmpi(ip.Results.Normalization, 'analytical')
    Av = Av ./ (2*pi*sv.^2);
end

switch ip.Results.Border
    case 'padded'
        for k = 1:np
            xw = -wv(k,1):wv(k,1);
            zw = -wv(k,2):wv(k,2);
            xa = xi(k)+xw;
            ya = yi(k)+xw;
            za = zi(k)+zw;
            [xg,yg,zg] = meshgrid(xw,xw,zw);
            g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*sv(k)^2)) .* exp(-(zg-dz(k)).^2 / (2*rv(k)^2));
            if strcmpi(ip.Results.Normalization, 'sum')
                g = Av(k)*g/sum(g(:));
            else
                g = Av(k)*g;
            end
            frame(ya,xa,za) = frame(ya,xa,za) + g;
        end
    case 'periodic'
        lbx = xi-wv(:,1);
        ubx = xi+wv(:,1);
        lby = yi-wv(:,1);
        uby = yi+wv(:,1);
        lbz = zi-wv(:,2);
        ubz = zi+wv(:,2);
        
        for k = 1:np
            shifts = [0 0 0];
            if lbx(k)<1
                shifts(2) = 1-lbx(k);
            elseif ubx(k)>nx
                shifts(2) = nx-ubx(k);
            end
            if lby(k)<1
                shifts(1) = 1-lby(k);
            elseif uby(k)>ny
                shifts(1) = ny-uby(k);
            end
            if lbz(k)<1
                shifts(3) = 1-lbz(k);
            elseif ubz(k)>nz
                shifts(3) = nz-ubz(k);
            end
            
            wi = -wv(k,1):wv(k,1);
            wiz = -wv(k,2):wv(k,2);
            [xg,yg,zg] = meshgrid(wi,wi,wiz);
            g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*sv(k)^2)) .* exp(-(zg-dz(k)).^2 / (2*rv(k)^2));
            if strcmpi(ip.Results.Normalization, 'sum')
                g = Av(k)*g/sum(g(:));
            else
                g = Av(k)*g;
            end
            xa = (xi(k)-wv(k,1):xi(k)+wv(k,1)) + shifts(2);
            ya = (yi(k)-wv(k,1):yi(k)+wv(k,1)) + shifts(1);
            za = (zi(k)-wv(k,2):zi(k)+wv(k,2)) + shifts(3);
            if all(shifts==0)
                frame(ya,xa,za) = frame(ya,xa,za) + g;
            else
                frame = circshift(frame, shifts);
                frame(ya,xa,za) = frame(ya,xa,za) + g;
                frame = circshift(frame, -shifts);
            end
        end
    case 'truncated'
        lbx = max(xi-wv(:,1),1);
        ubx = min(xi+wv(:,1),nx);
        lby = max(yi-wv(:,1),1);
        uby = min(yi+wv(:,1),ny);
        lbz = max(zi-wv(:,2),1);
        ubz = min(zi+wv(:,2),nz);
        
        for k = 1:np
            wx = (lbx(k):ubx(k)) - xi(k);
            wy = (lby(k):uby(k)) - yi(k);
            wz = (lbz(k):ubz(k)) - zi(k);
            [xg,yg,zg] = meshgrid(wx,wy,wz);
            g = exp(-((xg-dx(k)).^2+(yg-dy(k)).^2) / (2*sv(k)^2)) .* exp(-(zg-dz(k)).^2 / (2*rv(k)^2));
            if strcmpi(ip.Results.Normalization, 'sum')
                g = Av(k)*g/sum(g(:));
            else
                g = Av(k)*g;
            end
            xa = lbx(k):ubx(k);
            ya = lby(k):uby(k);
            za = lbz(k):ubz(k);
            frame(ya,xa,za) = frame(ya,xa,za) + g;
        end
end
