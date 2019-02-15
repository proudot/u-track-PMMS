function [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPosT,maximaAmpT,...
    bgAmpT,psfSigma,clusterPixels,firstFit)
%MMFINITGUESSLOWERUPPERBOUNDS calculates the initial guess and lower and upper bounds for mixture-model fitting
%
%SYNOPSIS [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPosT,maximaAmpT,...
%    bgAmpT,psfSigma,clusterPixels,firstFit)
%
%INPUT  maximaPosT : Particle positions.
%       maximaAmpT : Particle amplitudes.
%       bgAmpT     : Background amplitude.
%       psfSigma   : Gaussian sigma for approximating the point spread
%                    function.
%       clusterPixels: List of pixels belonging to cluster of local maxima.
%       firstFit   : 1 if this is the first time a group of local maxima is
%                    fitted, 0 otherwise.
%       PLEASE SEE detectSubResFeatures2D_V2 FOR PROPER CONTEXT.
%
%OUTPUT x0         : Initial guess.
%       lb         : Lower bound.
%       ub         : Upper bound.
%
%REMARKS This function in its current format is not really written for
%general use but as a function to be called by detectSubResFeatures2D_V2.
%
%Khuloud Jaqaman, August 2011

%feature positions
x0 = maximaPosT; %initial guess
lb = x0 - 2*psfSigma; %lower bound
minPos = min(clusterPixels);
lb(lb(:,1)<minPos(1),1) = minPos(1);
lb(lb(:,2)<minPos(2),2) = minPos(2);
if ~firstFit
    lb(end,:) = minPos;
end
ub = x0 + 2*psfSigma; %upper bound
maxPos = max(clusterPixels);
ub(ub(:,1)>maxPos(1),1) = maxPos(1);
ub(ub(:,2)>maxPos(2),2) = maxPos(2);
if ~firstFit
    ub(end,:) = maxPos;
end

%feature amplitudes
x0 = [x0 maximaAmpT];
lb(:,3) = eps;
ub(:,3) = 1;

%background intensity
x0 = x0';
x0 = [x0(:); bgAmpT];
lb = lb';
lb = [lb(:); eps];
ub = ub';
ub = [ub(:); 1];
