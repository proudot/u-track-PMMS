% Francois Aguet, 03/02/2012

function testADTestK()

% Example given in Table 4 of
% [1] Scholz and Stephens, J. Am. Stat. Assoc. 82(399), 1987

samples = {[38.7 41.5 43.8 44.5 45.5 46.0 47.7 58.0],...
           [39.2 39.3 39.7 41.4 41.8 42.9 43.3 45.8],...
           [34.0 35.0 39.0 40.0 43.0 43.0 44.0 45.0],...
           [34.0 34.8 34.8 35.4 37.2 37.8 41.2 42.8]};
[hval T cval] = adtestk(samples)
disp('-----------');


nexp = 1e4;
N = 1e2;
hvalAD = zeros(nexp,1);
for k = 1:nexp
    s1 = randn(1,N);
    s2 = randn(1,N);
    hvalAD(k) = adtestk({s1,s2});
end
sum(hvalAD)/nexp

disp('-----------');

% dx = 1;
% x = -10:dx:10;
% figure;
% hold on;
% for k = 1:numel(samples)
%     ni = hist(samples{k}, x);
%     ni = ni/sum(ni)/dx;
%     plot(x, ni, '.-', 'Color', rand(1,3));
% end


% N = 1e4;
% ns = 20;
% alpha = 0.05;
% 
% H_AD = NaN(1,N);
% H_KS = NaN(1,N);
% P_KS = NaN(1,N);
% H_L = NaN(1,N);
% 
% for k = 1:N
%     x = randn(1,ns);
%     x = (x-mean(x))/std(x);
%     
%     H_AD(k) = adtest(x, alpha);
%     [H_KS(k), P_KS(k)] = kstest(x, [], alpha);
%     H_L(k) = lillietest(x, alpha);
%     
% end
% 
% sum(H_AD)/N
% sum(H_KS)/N
% sum(H_L)/N

%%
% Exponential distributions

rng('default')
N = 1e4;
ns = 20;
hval = zeros(1,N);

for k = 1:N
    k
    hval(k) = adtest(-1*log(1-rand(1,ns)), 'Distribution', 'exponential');
end
sum(hval)/N


% test for a truncated distribution
hval = zeros(1,N);
mu = 1;
T = 0.2;
for k = 1:N
    samples = -mu*log(1-rand(1,ns));
    % truncate
    samples(samples<T) = [];
    samples = samples-T;
    hval(k) = adtest(samples, 'Distribution', 'exponential');
end
sum(hval)/N







