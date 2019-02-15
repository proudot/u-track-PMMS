function [sp,errStd,errV,corrS] = lsRandomFit(v1,v2,varargin)
%lsRandomFit: Spline (or linear) least-sqaure fitting between two random
%             variables with a large number of samples.
%
% Model: The goal of this function is to determine whether and how one random
%        variable depends on the other according to the following model:
%                         v2 = f(v1) + errV.
%
% SYNOPSIS:
%    [sp,errStd,errV] = lsRandomFit(v1,v2)
%    [sp,errStd,errV,corrS] = lsRandomFit(v1,v2,varargin)
%
% INPUT:
%    v1, v2: Vector of samples of two random variables.
%
%    Optional parameter/value pairs:
%       'model' : 'linear' or 'bspline'. When it is 'linear', a linear 
%                 least-square fitting is calculated. The output is a line.
%       'figH'  : If a figure handle is given, the fitting curve will be
%                 plotted.
%
% OUTPUT:
%    sp     : The sp-form of the fitting B-spline. If 'linear' is on, a line is 
%             given by sp = [a b] where v2 = a*v1+b+errV. Since we fit both
%              'v2' as a function of 'v1' and 'v1' as a function of 'v2', it is
%            really a cell array of sp-forms or lines.
%    errV   : A two-column matrix. The first column contains the error vector
%             when we fit 'v2' as a function of 'v1', the second when 'v1' 
%             as a function of 'v2'.
%    errStd : The standard deviation of the error vector 'errV'.
%    corrS  : Correlation score of the two random variables: 'v1' and
%             'v2'. When linear fitting is chosen, we give two
%             scores. One is the classic correlation coefficient. The
%             other is the cosine of the angle between the two
%             fitting lines.

if nargin < 2 | rem(nargin-2,2) ~= 0
   error('Wrong number of input arguments.');
end

model            = 'bspline';
figH             = [];

if nargin > 2 
   for k = 1:2:nargin-2
      switch varargin{k}
         case 'model'
            model = varargin{k+1};
         case 'figH'
            figH = varargin{k+1};
            if ~ishandle(figH)
               error('''figH'' is not a valid handle.');
            end
      end
   end
end

if ~strcmp(model,'linear') & ~strcmp(model,'bspline')
   error('The specified fitting model is not recognized.');
end

if length(v1) < 3 | length(v2) < 3
    sp     = [];
    errStd = [];
    errV   = [];
    corrS  = [];
    
    return;
end

[m,n] = size(v1);
if n == 1
   v1 = v1.';
end
[m,n] = size(v2);
if n == 1
   v2 = v2.';
end

%The minimum and maximum of 'v1' gives the range of fitting.
minV1 = min(v1);
maxV1 = max(v1);
minV2 = min(v2);
maxV2 = max(v2);

corrM = corrcoef(v1,v2);
corrS = corrM(1,2);

%%%%%% Calculate the fitting curve. We use least-square B-spline.
numKnots = 4;
errV     = NaN*ones(2,length(v2));
if strcmp(model,'bspline')
   spOrder  = 4;

   %Create the knot sequence.
   knots1 = augknt(linspace(minV1,maxV1,numKnots),spOrder);
   knots2 = augknt(linspace(minV2,maxV2,numKnots),spOrder);

   %Least sqaures spline fitting.
   sp1 = spap2(knots1,spOrder,v1, v2);
   sp2 = spap2(knots2,spOrder,v2, v1);

   errV(:,1) = v2 - fnval(sp1,v1);
   errV(:,2) = v1 - fnval(sp2,v2);
else
   %Linear least-square fitting.
   M = [v1.' ones(length(v1),1)];
   sp1 = M.'*M\(M.'*v2.');
   M = [v2.' ones(length(v2),1)];
   sp2 = M.'*M\(M.'*v1.');

   errV(1,:) = v2 - sp1(1)*v1 - sp1(2);
   errV(2,:) = v1 - sp2(1)*v2 - sp2(2);
   line1 = [1 sp1(1)];
   line2 = [sp2(1) 1];
   corrS(2) = sum(line1.*line2)/norm(line1)/norm(line2);
end

errStd(1) = std(errV(1,:));
errStd(2) = std(errV(2,:));
sp = {sp1,sp2};

if ishandle(figH)
   figure(figH); hold off;

   plot(v1,v2,'.'); hold on;

   plotellipse(eCenter,outlierThreshold*r1,outlierThreshold*r2,eAngle,'m');

   %Plot the fitting curve.
   plotRange1 = linspace(minV1,maxV1,5*numKnots);
   if strcmp(model,'bspline')
      plot(plotRange1,fnval(sp1,plotRange1),'r','lineWidth',3);
   else
      plot(plotRange1,sp1(1)*plotRange1+sp1(2),'r','lineWidth',3);
   end
   
   %Plot the fitting curve.
   plotRange2 = linspace(minV2,maxV2,5*numKnots);
   if strcmp(model,'bspline')
      plot(fnval(sp2,plotRange2),plotRange2,'g','lineWidth',3);
   else
      plot(sp2(1)*plotRange2+sp2(2),plotRange2,'g','lineWidth',3);
   end
end


