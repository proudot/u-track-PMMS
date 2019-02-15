%adtestkiterative(samples) iteratively computes the k-sample Anderson-Darling test on a set of samples.
% The test is computed on iteratively smaller subsets, until a subset of samples that
% pass the null hypothesis (samples are from the same distribution) is found.
%
% Inputs:
%    samples : cell array of samples
%
% Outputs:
%         id : index of samples, empty if no combination passes the test

% Francois Aguet, 05/24/2012

function id = adtestkiterative(samples)
id = [];
ns = numel(samples);

for n = ns:-1:2
    combs = nchoosek(1:ns,n);
    nc = size(combs,1);
    for c = 1:nc
        h = adtestk(samples(combs(c,:)));
        if h==0
            id = combs(c,:);
            return
        end
    end
end
