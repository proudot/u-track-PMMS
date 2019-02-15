function [listOfFiles,directory] = searchFiles2(includeString1,includeString2,excludeString,directory,includeSubDirectories,selectionMode,dirStr)
%searchFiles2 is an utility to search for files containing 1-2 specific strings
%
%SYNOPSIS [listOfFiles,directory] = searchFiles(includeString1,includeString2,excludeString,directory,includeSubDirectories,selectionMode)
%
%INPUT    includeString1: string contained in the filenames you are looking for
%         includeString2 (opt): another string contained in the filenames you are looking for  
%         excludeString (opt): string not contained in the filenames you are looking for
%         directory (opt): directory to search. if empty, current directory is searched (default)
%                          if 'ask', program asks for directory
%         includeSubDirectories (opt): whether to search subdirectories or not (0/{1})
%         selectionMode (opt): which file(s) to select if there are several files matching the
%                              includeString within one directory
%                              {'all'}-all; 'new'-newest; 'old'-oldest; 'GUI'-open GUI to select one file
%         dirStr (opt): string for directory search box (default 'select a
%                       directory to search')
%
%OUTPUT  listOfFiles: cell array with {[fileName] [directoryName]}
%        directory: path of directory searched
%
%c: 9-05 kathryn - added possibility to look for two strings in file name and specify dirStr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---test input---
if nargin<1|isempty(includeString1)
    error('not enough input arguments or empty includeString1')
end

%includeString1
if ~isstr(includeString1)
    error('includeString1 has to be a string')
end

%includeString2
if nargin>1&~isempty(includeString2)
    if ~isstr(includeString2)
        error('includeString2 has to be a string')
    end
else
    includeString2 = [];
end

%excludeString
if nargin>2&~isempty(excludeString)
    if ~isstr(excludeString)
        error('excludeString has to be a string')
    end
else
    excludeString = [];
end

%dirStr
if nargin<7|isempty(dirStr)
    dirStr = 'select a directory to search'
end

%directory (ask if necessary)
if nargin>3
    if isempty(directory)
        directory = pwd;
    elseif strcmp(directory,'ask')
        directory = uigetdir(pwd,dirStr);
        if directory == 0
            error('searchFiles aborted by user')
        end
    elseif ~isdir(directory)
        error([directory,' is not a valid directory!'])
    end
end
if nargin<4|isempty(directory)
    directory = pwd;
end

%includeSubDirectories
if nargin<5|isempty(includeSubDirectories)
    includeSubDirectories = 1;
end

%selectionMode
if nargin<6|isempty(selectionMode)
    selectionMode = 'all';
else
    if ~(strcmp(selectionMode,'new')+strcmp(selectionMode,'old')+strcmp(selectionMode,'GUI')+strcmp(selectionMode,'all'))
        error('wrong selectionMode')
    end
end

%---end test input---


%---collect files---

%look for directories and wanted files in current directory, then check all subdirs

%init variables

dirs2checkInitLength = 100;
dirs2checkLength     = 100;
dirs2check           = cell(100,1);

dirs2check{1} = directory;
dirs2checkCt  = 1;
topDir2check  = 1;

listOfFilesInitLength = 100;
listOfFilesLength     = 100;
listOfFiles           = cell(100,2);
listOfFilesCt         = 0;


while ~isempty(dirs2check) && ~isempty(dirs2check{topDir2check})
    %init/empty var
    pathCell = {};
    
    %read currentDir
    currentDir = dirs2check{topDir2check};
    
    %list all files of current directory
    currentDirList = dir(currentDir);
    
    %look for subdirectories and add to dirs2check
    isDirList = cat(1,currentDirList.isdir);
    if length(isDirList)>2
        subDirIdx = find(isDirList(3:end))+2;
        for i=1:length(subDirIdx)
            
            % we will add a directory - count how many this will make
            dirs2checkCt = dirs2checkCt + 1;
            
            % make sure that the dirList is long enough - make sure we do
            % not get problems with topDir2check
                if dirs2checkCt + 1 > dirs2checkLength
                    tmpDirs2check = dirs2check;
                    newDirs2checkLength = dirs2checkLength + dirs2checkInitLength;
                    dirs2check = cell(newDirs2checkLength,1);
                    dirs2check(1:dirs2checkLength) = tmpDirs2check;
                    dirs2checkLength = newDirs2checkLength;
                end
            
            dirs2check{dirs2checkCt} = [currentDir,filesep,currentDirList(subDirIdx(i)).name];
        end
    end %if length(isDirList)>2
    
    %look for files in current directory and store them
    if isstr(includeString2)
        newFiles = intersect(chooseFile(includeString1,currentDir,selectionMode,excludeString),...
            chooseFile(includeString2,currentDir,selectionMode,excludeString));
    else
        newFiles = chooseFile(includeString1,currentDir,selectionMode,excludeString);
    end
        
    if ~isempty(newFiles)
        if ~iscell(newFiles)
            newFiles = cellstr(newFiles);
        end
        [pathCell{1:size(newFiles,1)}] = deal(currentDir);
        
        % we will add a file - count how many this will make
        numNewFiles = length(newFiles);
        listOfFilesCtStart = listOfFilesCt + 1;
        listOfFilesCt = listOfFilesCt + numNewFiles;
        
        % make sure that the dirList is long enough
        if listOfFilesCt > listOfFilesLength
            tmpListOfFiles = listOfFiles;
            newListOfFilesLength = listOfFilesLength + listOfFilesInitLength;
            listOfFiles = cell(newListOfFilesLength,2);
            listOfFiles(1:listOfFilesLength,:) = tmpListOfFiles;
            listOfFilesLength = newListOfFilesLength;
        end
        
        listOfFiles(listOfFilesCtStart:listOfFilesCt,:) = [newFiles, pathCell'];
        
    end %if ~isempty(newFiles)
    
    % move on to next directory
    topDir2check = topDir2check + 1;
    
    %check wheter we want to look at subDirs
    if ~includeSubDirectories
        dirs2check = {};
    end
    
    
end %while ~isempty(dirs2check)

%---end collect files---

% remove placeholders
listOfFiles(listOfFilesCt+1:end,:)=[];