function filename = enumFileName(directory, filename, type)
% Find all .mat files in input directory and return a new file name with 
% proper enumeration (e.g filename1, filename2, ...)
%
% Input: 
%   directory: the directory for the new file
%   filename: the base of file name
%   type: file type, e.g. .mat, .mov, etc
%
% Output: 
%   filename: full name. E.g. filename1, filename2, ...
%
%
% Chuangang Ren
% 12/2010

if nargin < 3
    type = '.mat';
end

if ~exist(directory,'dir')
    inputbase = regexp(filename, type,'split');
    inputbase = inputbase{1};
    filename = [inputbase '1' type];
    return
end

if strcmp(directory(end), filesep)
    matFileName = dir([directory '*' type]);
else
    matFileName = dir([directory filesep '*' type]);
end
matFileName = arrayfun(@(x)(x.name),matFileName,'UniformOutput',false);

[x1 base digits4Enum x4] = cellfun(@(x)getFilenameBody(x), matFileName,'UniformOutput',false);

inputbase = regexp(filename, type,'split');
inputbase = inputbase{1};

new = 0;

if isempty(base)
    
    new = 1;
   
else
   
    samebase = strcmp(inputbase, base);
    if ~any(samebase)
        
        new = 1;
        
    else
       
        dig = max(str2double(digits4Enum(samebase)));
        if isnan(dig)
            
            new = 1;
        else
            filename = [inputbase num2str(dig+1) type];
        end
        
    end
        
end

if new
    filename = [inputbase '1' type];
end
