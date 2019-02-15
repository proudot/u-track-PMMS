function [mdsX,mdsY] = modeFinder(map,type,parameter)
%MODEFINDER seeks for statistically significant modes in a 2D arry (e.g. a cross correlation map)
% 
% SYNOPSIS [mdsX,mdsY] = modeFinder(map,type,parameter)
%
% INPUT map  : 2D array with data
%       type : either be 'fixed', for which the parameter specififies how many mode are extracted
%                     or 'automatic', for which the software tries to retrieve any statistically
%                                     significant mode
%                     or 'fixedButCheck'  (description see below)
%
%              For the type 'automatic' the parameter represents the threshold testing the following ratio:
%
%              max difference in a sorted locmax vector
%              ----------------------------------------------------------
%              standard deviation of all the nonsignificant locmax values
%
%              This tests the statistical significance of the modes that seem do dominate the map
%              For this calculation mode two parameters are required
%
%              parameter(1) : ratio threshold
%              parameter(2) : median filter size, use to remove a low frequency component in the mode map
%
%                       'fixedButCheck' combines the two before mentioned types. However, if more than the 
%                       specified number of modes are detected as significant, then the 
%                       set is reduced to the number of desired modes. If the opposite is the
%                       case, i.e. less modes than the desired number are significant,
%                       then only the significant modes are returned
%                       For this mode we need 3 parameters;
%                       parameter(1) : number of sought modes
%                       parameter(2) : ratio threshold
%                       parameter(3) : median filter size (see above)
        

calcType = type;

switch type
case 'fixed',
    nModes = parameter(1);
case 'automatic',
    if(length(parameter) ~= 2)
        error('for the FIXEDBUTCHECK option, 2 parameters have to be specified');
    end;
    thresh = parameter(1);
    medFiltSize = parameter(2);
case 'fixedButCheck',
    calcType = 'automatic';
    if(length(parameter) ~= 3)
        error('for the FIXEDBUTCHECK option, 2 parameters have to be specified');
    end;
    nModes = parameter(1);
    thresh = parameter(2);
    medFiltSize = parameter(3);
end


if strcmp(calcType,'automatic')
    % there might be a large low frequency component in the map. To perform 
    % statistics on the variation of the local maxima peak we have to remove this
    % from the data
    map = map - medfilt2(map,[medFiltSize medFiltSize], 'symmetric');    
end;

% calculate locmax values
lMaxMap = locmax2d(map,[3 3]);

% retrieve all locmax values in a list and sort them
[posRow,posCol,vals]=find(lMaxMap);
[srtVals,iSrtVals] = sort(vals);

switch calcType
case 'fixed', 
    mdsX = posCol(iSrtVals(end:-1:end-nModes+1));
    mdsY = posRow(iSrtVals(end:-1:end-nModes+1));
case 'automatic',
    % find the maximum step in the sorted list; values beyond
    % that index may be significant modes; values below that index must
    % be considered as arbitrary random locmax values
    dSrtVals = diff(srtVals);
    [dummy,iMaxDiff] = max(dSrtVals);
    
    % the step size has to be significantly larger than the variation of 
    % the random locmax values; otherwise, there is no meaningful mode that clearly
    % stands out of the rest of the values
    if((dSrtVals(iMaxDiff)/std(srtVals(1:iMaxDiff))) > thresh)
        mdsX = posCol(iSrtVals(length(srtVals):-1:iMaxDiff+1));
        mdsY = posRow(iSrtVals(length(srtVals):-1:iMaxDiff+1));
    else
        mdsX=[];
        mdsY=[];
    end;

end;

if(strcmp(type,'fixedButCheck'))
    if(length(mdsX) > nModes)
       mdsX = mdsX(1:nModes);
       mdsY = mdsY(1:nModes);
   end;
end;


% subpixel positioning of the modes
[peakDX, peakDY] = meshgrid(-1:1,-1:1);
for(nM = 1: length(mdsX))
    peakD = map(mdsY(nM)-1:mdsY(nM)+1,mdsX(nM)-1:mdsX(nM)+1);
    % biquadratic fit to the peak
    a = biquadfit(peakDX,peakDY,peakD);
    % equation system setting up the peak finder based on the biquadratic 
    % polynom
    A = [2*a(4),a(6);a(6),2*a(5)];
    b = -[a(2);a(3)];
    shift = A\b;
    mdsX(nM)= mdsX(nM) + shift(1);
    mdsY(nM)= mdsY(nM) + shift(2);
end;