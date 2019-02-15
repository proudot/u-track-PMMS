function mkClrDir(dirPath)
%MKCLRDIR makes sure that the specified directory exists AND is empty 
% 
% This is just a little function for creating / settin up output
% directories. It checks if a directory exists and if not, makes it. If the
% directories does exist and contains files, these are deleted.
% 
% Input:
% 
%   dirPath - the path to the directory to make/clear.
% 
% Hunter Elliott
% 6/2010

if nargin < 1 || isempty(dirPath)
    error('You must specify the directory to set up!')
end

if ~exist(dirPath,'dir')
    mkdir(dirPath)
else
    %Check for files in the directory
    inDir = dir([dirPath filesep '*']);
    if ~isempty(inDir)
        %Remove the . and .. from dir (on linux)
        inDir = inDir(arrayfun(@(x)(~strcmp('.',x.name) ...
            && ~strcmp('..',x.name)),inDir));
        arrayfun(@(x)(delete([dirPath filesep x.name])),inDir);
    end
end