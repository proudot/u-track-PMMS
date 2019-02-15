function [filePath, fileName] = optionalFileInput(fileIn,filterString,promptString)
%OPTIONALFILEINPUT handles optional file name input or user prompting with dialogue
%
% [filePath, fileName] = optionalFileInput;
% [filePath, fileName] = optionalFileInput(fileNameIn);
% [filePath, fileName] = optionalFileInput([],filterString,promptString)
% [filePath, fileName] = optionalFileInput(fileNameIn,filterString,promptString)
%
%   Just a simple function to handle simple tasks related to having a file
%   name either be input directly into a function or if not prompting the
%   user to select a file with a file selection GUI dialogue box.
%
% Hunter Elliott
% 8/2013

if nargin < 2
    filterString = {'*.*'};
elseif ~iscell(filterString)
    filterString = {filterString};
end

if nargin < 3 || isempty(promptString)
    promptString = 'Please select a file:';
end

if nargin < 1 || isempty(fileIn)    
    [fileName,filePath] = uigetfile(filterString(:),promptString);
    if fileName == 0
        return
    end
else    
    [filePath,fileName,fExt]= fileparts(fileIn);        
    fileName= [fileName fExt];
    if isempty(filePath)
        filePath = pwd;
    end
    filePath = [filePath filesep];
end