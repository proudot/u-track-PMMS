%[Y] = WSTARCK2(IMG,K,IRECON) denoises the input image  by hard thresholding in the wavelet subbands.
%
% Inputs:
%       img : input image   
%         k : decomposition level of the wavelet transform
%    irecon : flag to include (1) or not (0) the low-pass component in the reconstruction
%
% Reference:  
% [1] Olivo-Marin, "Extraction of spots in biological images using multiscale products," Pattern Recoginition 35, pp. 1989-1996, 2002.
%
% Note: The implementation is incorrect, the thresholding mask should be calculated
%       separately for each scale instead of being cumulative. Since this code has
%       been used for published work, and the bug has effects on downstream functions,
%       it is left as is. -FA

% Shann-Ching Sam Chen, 2008. Last modified by Francois Aguet.

function y = wstarck2(img, k, irecon)

[ny,nx] = size(img);

y = zeros(ny,nx);
s = zeros(ny,nx);

% Calculate 'a trous' wavelet transform (spline wavelets)
W = awt(img, k);

for ind = 1:k
    tmp = W(:,:,ind);
    s(abs(tmp)>=3*std(tmp(:))) = 1;
    y = y + s.*tmp;
end

if irecon==1
    y = y + W(:,:,end);
end
