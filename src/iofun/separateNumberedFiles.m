function separateNumberedFiles(parentDir,numString,nameString,dryRun)
%SEPARATENUMBEREDFILES separates numbered files into their own numbered directories
% 
% separateNumberedFiles(parentDir,numString)
% separateNumberedFiles(parentDir,numString,nameString,testRun)
% OK = separateNumberedFiles(...);
%
% This function is designed to take files which are numbered somewhere in
% their file name and put them into their own, numbered directories. This
% is generally useful for raw microscopy data, where several
% movies/experiments may be stored in the same directory but need to be
% separated. These files often have more than one number in their name,
% which necessitates specifying the location of the numbering of interest.
%
% For example, if a folder named "exampleDir" contains the files:
% 
%   movie1_t1.tif
%   movie1_t2.tif
%   movie1_t3.tif
%   movie2_t1.tif
%   movie2_t2.tif
%   movie2_t3.tif
% 
% We want to separate based on the movie number, not the timepoint number.
% Therefore, the command
% 
%   separateNumberedFiles(parentDir,'movie#')
% 
% Would put the first 3 files in a folder called "movie1" and the last
% three files in a folder called "movie2", while the command
% 
%   separateNumberedFiles(parentDir,'t#')
%
% would put the first and 4th files in one directory named "t1", the second
% and fifth in a directory named "t2" and the 3rd and 6th in a directory
% named "t3"
%
% Input:
%
%   parentDir - The directory containing the files to separate. Optional.
%   The default is to use the current working directory.
%
%   numString - A character array specifying the location of the numbering
%   within the file names with the # character. See the example above for
%   more explanation. If each file only has one number in its name, this is
%   not required. If the string matches more than one number in the file
%   name, the first number will be used.
%
%   nameString - A character array specifying the names to use for the
%   folders, also containing a # character to designate where the number
%   goes. Optional. If not input, the folders will be named based on the
%   numString input.
%
%   dryRun - Use to test prior to actually moving files. If true, no files
%   will be copied, but the list of what would have been moved will still
%   be displayed.
%   Optional. Default is false.
%
% Hunter Elliott
% 11/2010
%

%% ------- Parameters ------- %%

ns ='#';%String used to designate numberlocations

%% -------- Input --------- %%

if nargin < 1 || isempty(parentDir)
    parentDir = pwd;
elseif ~exist(parentDir,'dir')
    error('Specified parent directory is not a valid directory!')
end
    
if nargin < 2 || isempty(numString)
    numString = ns;
elseif isempty(regexp(numString,ns,'ONCE')) || numel(regexp(numString,ns)) > 1
    error('The numString input must contain exactly one # character!')
end

if nargin < 3 || isempty(nameString)
    nameString = numString;
elseif isempty(regexp(nameString,ns,'ONCE')) || numel(regexp(nameString,ns)) > 1
    error('The nameString input must contain exactly one # character!')
end

if nargin < 4 || isempty(dryRun)
    dryRun = false;
end

%% --------- Init --------- %%

%Find all files in the directory
fileNames = dir(parentDir);
%Remove directories from the list
fileNames = fileNames(arrayfun(@(x)(~x.isdir),fileNames));

if isempty(fileNames)
    error('No files were found in the specified parent directory!')
end

%Get index of the # in numString
iPs = regexp(numString,ns);

%Convert numString to regular expression
numString = regexprep(numString,ns,'\\d+');

nFiles = numel(fileNames);

%% ------ Separate Files ----- %%


%First we need to determine which numbers are present and extract them from
%the file names

fileNumbers = nan(nFiles,1);

for j = 1:nFiles

    %Find the first place the numbering string matches in this file name
    iNS = min(regexp(fileNames(j).name,numString));
    
    if ~isempty(iNS)
        %Find the first non-numeric element after this to determine the number
        %of digits        
        nDig = min(regexp(fileNames(j).name(iNS+iPs-1:end),'\D'))-1;            
        
        %Convert the number string to a number
        fileNumbers(j) = str2double(fileNames(j).name(iNS+iPs-1:iNS+iPs-1+nDig-1));
    
        
    end        

end

%Make sure at least one file matched the numbering string
if ~any(~isnan(fileNumbers))
    error('No files in the specified directory matched the specified numString!')
elseif numel(unique(fileNumbers(~isnan(fileNumbers)))) == 1
    error('All of the files matching the numbering string have the same number!')
end

%Get the numbers for the directories to be used
dirNums = unique(fileNumbers(~isnan(fileNumbers)));
nDir = numel(dirNums);

disp(['Found ' num2str(nDir) ' different groups of numbered files.'])

%Go through these directories, create them and put the files in them
for j = 1:nDir
    
    %Create this directory
    dirName = regexprep(nameString,'#',num2str(dirNums(j)));    
    
    if ~dryRun
        mkdir([parentDir filesep dirName])
    end
    
    iThisDir = find(fileNumbers == dirNums(j));
    nFiles = numel(iThisDir);
    
    if ~dryRun
        
        disp(['Moving ' num2str(nFiles) ' files into directory "' dirName '"'])

        for k = 1:nFiles

            movefile([parentDir filesep fileNames(iThisDir(k)).name],...
                 [parentDir filesep dirName filesep fileNames(iThisDir(k)).name]);


        end 
    else
        disp(['Dry run, not doing anything. Would have moved ' num2str(nFiles) ' files into directory "' dirName '"'])        
    end
end



