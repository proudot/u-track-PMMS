function [framesNoCands,errFlag] = checkCandsStatus(candsParam)
%CHECKCANDSSTATUS find frames where all cands have status 0 (i.e. are insignificant)
%
%SYNOPSIS [framesNoCands,errFlag] = checkCandsStatus(candsParam)
%
%INPUT  candsParam    : Structure with fields
%           .candsDir     : Directory where cands (initial maxima) are stored.
%           .firstCandsNum: Numerical index of first cands file.
%           .lastCandsNum : Numerical index of last cands file.
%           .digits4Enum  : Number of digits used to enumerate cands files.
%
%OUTPUT framesNoCands : Array indicating frames where all cands have status
%                       0 (i.e. empty frames).
%
%Khuloud Jaqaman, December 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

framesNoCands = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--checkCandsStatus: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cands files parameters
candsDir = candsParam.candsDir;
firstCandsNum = candsParam.firstCandsNum;
lastCandsNum = candsParam.lastCandsNum;
digits4Enum = candsParam.digits4Enum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Status checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assign leading zeros in numerical index of cands files
switch digits4Enum
    case 4
        leadingZeros(1).value = '000';
        leadingZeros(2).value = '00';
        leadingZeros(3).value = '0';
        leadingZeros(4).value = '';
    case 3
        leadingZeros(1).value = '00';
        leadingZeros(2).value = '0';
        leadingZeros(3).value = '';
    case 2
        leadingZeros(1).value = '0';
        leadingZeros(2).value = '';
    case 1
        leadingZeros(1).value = '';
end

%go through the cands files

for i=min(9999,lastCandsNum):-1:max(1000,firstCandsNum)

    %load cands
    eval(['load ' candsDir 'cands' leadingZeros(4).value num2str(i) ';'])

    %get cands status
    status = vertcat(cands.status);

    %find cands with status 1
    indx = find(status==1);

    %if there aren't any, report frame
    if isempty(indx)
        framesNoCands = [framesNoCands; i];
    end

end

for i=min(999,lastCandsNum):-1:max(100,firstCandsNum)

    %load cands
    eval(['load ' candsDir 'cands' leadingZeros(3).value num2str(i) ';'])

    %get cands status
    status = vertcat(cands.status);

    %find cands with status 1
    indx = find(status==1);

    %if there aren't any, report frame
    if isempty(indx)
        framesNoCands = [framesNoCands; i];
    end

end

for i=min(99,lastCandsNum):-1:max(10,firstCandsNum)

    %load cands
    eval(['load ' candsDir 'cands' leadingZeros(2).value num2str(i) ';'])

    %get cands status
    status = vertcat(cands.status);

    %find cands with status 1
    indx = find(status==1);

    %if there aren't any, report frame
    if isempty(indx)
        framesNoCands = [framesNoCands; i];
    end

end

for i=min(9,lastCandsNum):-1:max(1,firstCandsNum)

    %load cands
    eval(['load ' candsDir 'cands' leadingZeros(1).value num2str(i) ';'])

    %get cands status
    status = vertcat(cands.status);

    %find cands with status 1
    indx = find(status==1);

    %if there aren't any, report frame
    if isempty(indx)
        framesNoCands = [framesNoCands; i];
    end

end


%%%%% ~~ the end ~~ %%%%%
