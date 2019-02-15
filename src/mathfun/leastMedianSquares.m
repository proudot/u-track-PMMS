function [u,inlierIdx,sigma0,sigmaU,inlierResiduals] = leastMedianSquares(...
    functionString,u0,options,parameters)
%leastMedianSquare calculates the least median squares for a function handed down as string using parameters (including x/y-data) from the structure parameters
%
%SYNOPSIS [u,inlierIdx,sigma0,sigmaU,inlierResiduals] = leastMedianSquares(...
%    functionString,u0,options,parameters)
%
%INPUT functionString       string which specifies the function to be
%                           minimized, e.g. '(y-(u(1)*x+u(2)))'. The
%                           unknown has to be 'u'. Do not square!
%      u0                   vector of the initial values for the unknowns in the function
%      options              option structure for optimization routines. Set them with optimset
%      parameters           structure containing all parameters of the
%                           function e.g. y and x as parameters.x and
%                           parameters.y. Note that leastMedianSquares will
%                           attempt to remove outlier-rows in all
%                           parameters with length > 1. To avoid that you can
%                           use inlierIdx and nInlierIdx as parameters
%                           inside your function.
%                           Initialize nInlierIdx as nDatapoints, and
%                           inlierIdx as 1:nDatapoints.
%
%AS TESTING THE INPUT FOR SPELLING ERRORS IN THE FUNCTION STRING OR THE FIELDNAMES OF PARAMETERS IS FAR TOO COMPLICATED
%TO DO WITHIN THE FUNCTION, PLEASE BE CAREFUL WITH YOUR INPUT!
%
%OUTPUT  u                  vector with estimates for the unknowns u, e.g. a and b
%        inlierIdx          rows in the data that are inliers (get
%                           outlierIdx with setdiff(1:nData,inlierIdx) )
%        sigma0             estimate for std of error of the data according to the fit
%        sigmaU             vector with std of the parameters according to
%                           lsqnonlin
%        inlierResiduals    vector of residuals from the fit, naturally
%                           only for inlier points
%
%REMEMBER TO CHECK THE NUMBER OF LOCAL MINIMA IN THE REGION OF YOUR
%SOLUTION (lsqnonlin is sensitive to intial conditions)
%
% Possible improvements: - Use weights (e.g. field in parameter structure),
%                          and build the objective function anlalogous to
%                          linearLeastMedianSquares
%                       -  Allow boundaries (switch to fminbnd in that
%                          case)
%
%c: Jonas, 3/04
% 11/07 made code robust to NaNs - also: use lsqnonlin. That makes the code
% a little less robust, but it returns much better solutions and allows to
% return uncertainties.
%   10/09: Code now outputs the residuals for inlier points - KJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note: this is an old function, but I won't rework it for prettiness

%initialize
parmString=[];
k=3; %value important for calculation of sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

u0 = u0(:); %#ok<NASGU>

%get names of parameters for fminsearch
parameterNames=fieldnames(parameters);
parameterNamesN=length(parameterNames);

%prepare fminsearch
funStringFmin = ['(' functionString,').*(',functionString ')']; % we want median(fs.*fs), not median fs'*fs!!
funStringMedian=['''nanmedian(',funStringFmin,')'', ''u'''];


for i=1:parameterNamesN
    parmString=[parmString,', ',parameterNames{i}]; %#ok<AGROW> % write variable into parameter string
    eval([parameterNames{i},'=parameters.',parameterNames{i},'(:);']); % assign value to variable and make col-vectors
    funStringMedian=[funStringMedian,', ''',parameterNames{i},'''']; %#ok<AGROW> % write variable into inline object
end

% fminsearch
u=eval(['fminsearch(inline(',funStringMedian,'),u0,options ',parmString,');']); %#ok<NASGU>
% debug
% uMedian = u;

%calculate statistics
eval(['res2=', funStringFmin,';']);
medRes2 = nanmedian(res2);

%testvalue to calculate weights
testValue=res2/(magicNumber2*medRes2);

%inlierIdx: weight 1, badRows: weight 0
% this will also remove any row that evaluates to NaN!
inlierIdx=find(testValue<=k^2);

%ssq=sum(res2);
sigma0=sqrt(sum(res2(inlierIdx))/(length(inlierIdx)-4));


% write function string for least squares formulation

%Jonas's original version:
% funStringLsq = ['''((' funStringFmin '))'', ''u'' ']; %now we want the sum of funStringFmin

%KJ: fixed bug in original version, where the squared function was passed
%to lsqnonlin, instead of the function itself.

%CORRECT VERSION:
funStringLsq = ['''((' functionString '))'', ''u'' '];

% update parameters: take only good rows - or write inlierIdx into
% parameters
inlierList = ~cellfun('isempty',regexpi(parameterNames,'inlier'));
isInlier = any(inlierList);
if isInlier
    parameters.inlierIdx = inlierIdx;
    parameters.nInlierIdx = length(inlierIdx);
end
for pn=parameterNames'
    pn=char(pn);
    % assign parameters

    lpn = length(parameters.(pn));
    if ~isInlier && lpn>1
        % assign NaN
        eval([pn,'=parameters.',pn,'(inlierIdx);']); % reassign
    else
        eval([pn,'=parameters.',pn,';']); % don't shorten
    end
    funStringLsq=[funStringLsq,', ''',pn,'''']; %#ok<AGROW> % write variable into inline object
end


% use uMedian to init fminsearch
% and replace fminsearch with lsqnonlin
[u,resnorm,residual,exitflag,output,lambda,jacobian]=...
    eval(['lsqnonlin(inline(',funStringLsq,'),u,[],[],options',parmString,');']);

% for dof: count also outliers
dof = length(res2)-length(u);
sigmaU = full(sqrt(diag(inv(jacobian'*jacobian))*resnorm/dof));

%save inlierResiduals for output
inlierResiduals = residual;
