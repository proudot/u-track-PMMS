function copySubsetOfFiles(dirSourceTarget,indices2copy,filenameBase,...
    filenameExt,digits4Enum,renumber)
%COPYSUBSETOFFILES copies a subset of files from one directory to another, with the possibility of renumbering
%
%SYNPOSIS copySubsetOfFiles(dirSourceTarget,indices2copy,filenameBase,...
%    filenameExt,digits4Enum,renumber)
%
%INPUT  dirSourceTarget: 2-by-1 cell array indicating source directory to
%                     copy files from and target directory to copy files
%                     to.
%       indices2copy: Row vector of file indices to copy.
%       filenameBase: File name base (everything before the numbers).
%       filenameExt : File name extension (everything after the ".").
%       digits4Enum : Number of digits used to enumerate files.
%       renumber    : 1 to renumber the copied files, starting with 1 and
%                     using increments of 1; 0 otherwise. When files are
%                     renumbered, "new_" is inserted after the filenameBase
%                     and before the numbers.
%                     Optional. Default: 0.
%
%OUTPUT NONE
%
%Khuloud Jaqaman, April 2010

%% Input

%check whether correct number of input arguments was used
if nargin < 5
    disp('--copySubsetOfFiles: Incorrect number of input arguments!');
    return
end

%get source and target directories
sourceDir = dirSourceTarget{1};
targetDir = dirSourceTarget{2};

%check whether to renumber files
if nargin < 6 || isempty(renumber)
    renumber = 0;
end

%% Action!

%to keep the user informed of progress...
progressText(0,'Copying files');
numFiles2copy = length(indices2copy);
iDummy = 0;

%based on renumbering choice ...
if renumber %if file names are renumbered
    
    %initialize index for renumbering
    iNew = 0;
    
    %go over files to be copied
    for iFile = indices2copy
        
        %get file name
        file2copy = [filenameBase sprintf(['%0' num2str(digits4Enum) 'i'],iFile) '.' filenameExt];
        
        %construct new file name
        iNew = iNew + 1;
        newName = [filenameBase 'new_' sprintf(['%0' num2str(digits4Enum) 'i'],iNew) '.' filenameExt];
        
        %copy file
        copyfile(fullfile(sourceDir,file2copy),fullfile(targetDir,newName));
        
        %to keep the user informed of progress...
        iDummy = iDummy + 1;
        progressText(iDummy/numFiles2copy,'Copying files');
        
    end
    
else %if file names stay the same
    
    %go over files to be copied
    for iFile = indices2copy
        
        %get file name
        file2copy = [filenameBase sprintf(['%0' num2str(digits4Enum) 'i'],iFile) '.' filenameExt];
        
        %copy file
        copyfile(fullfile(sourceDir,file2copy),targetDir);
        
        %to keep the user informed of progress...
        iDummy = iDummy + 1;
        progressText(iDummy/numFiles2copy,'Copying files');
        
    end
    
end

%% ~~~ the end ~~~
