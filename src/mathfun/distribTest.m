function [confValue]=distribTest(pop1,pop2)
% distribTest returns the percent confidence that two distributions are different using K-S test
%
% here we test whether two populations, pop1 and pop2, are different by
% bootstrapping. first pop1 is used to find a calibrated p-value: we sample
% n1 values (with replacement) twice, and compare the mean-subtracted
% sample distributions using a K-S test.  by doing this many times, we get
% a distribution of p-values. the 5th percentile represents the
% bootstrapped p-value, thresh.  then we perform the K-S test many times
% with pop1 and pop2.  the proportion of resulting p-values smaller than
% thresh represent the confidence that the two populations are in fact
% different.


nReps=1000;

n1=length(pop1);
n2=length(pop2);

p1=zeros(nReps,1);       
p2=zeros(nReps,1);
for i=1:nReps
    % sample both populations with replacement, but sample the first twice
    s1a=randsample(pop1,n1,true);
    s1b=randsample(pop1,n1,true);
    
    s2a=randsample(pop2,n2,true);
    
    [h p1(i)]=kstest2(s1a-mean(s1a),s1b-mean(s1b)); % calibration
    [h p2(i)]=kstest2(s1a-mean(s1a),s2a-mean(s2a));
end

thresh=prctile(p1,5);
confValue=100*(sum(p2<thresh)/length(p2));