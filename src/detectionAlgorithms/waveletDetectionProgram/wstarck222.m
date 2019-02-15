%[wm] = wstarck222(img, k) calculates the multiscale wavelet product.

% Shann-Ching Sam Chen, 2008 (?). Last modified by Francois Aguet.

function [wm] = wstarck222(img, k)

[ny,nx] = size(img);

% Calculate 'a trous' wavelet transform (spline wavelets)
W = awt(img, k);

wm = ones(ny,nx);
for ind = 1:k 
    wm = abs(wm.*W(:,:,ind));
end
