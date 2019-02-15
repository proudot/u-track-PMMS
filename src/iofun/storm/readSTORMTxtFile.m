function stormData = readSTORMTxtFile(fileName,varargin)
%READSTORMTXTFILE reads localization data from storm txt files of various formats
%
% stormData = readSTORMTxtFile
% stormData = readSTORMTxtFile(fileName)
%
%   This is a generic wrapper function for various storm txt file readers.
%   Currently supported formats are:
%
%       Nikon, Zeiss
%
%   These proprietary formats are read with dedicated functions and then
%   converted into a standardized format
%
% Output:
%
%   stormData - A standardized structure 

%Hunter Elliott
% 1 / 2015

%% --------- Parameters ------- %%

%I know, I should just create a storm data object but I'm in a hurry....
readerFun = {@readNikonSTORM,@readZeissSTORM};%Functions to read the files
testFun = {@isFileNikonSTORM,@isFileZeissSTORM};%Functions to check the format (since the file extensions are the same)
readerName = {'Nikon','Zeiss'};
nReader = numel(readerFun);

%% --------- Input --------- %%

if nargin < 1
    fileName = '';
end

[filePath,fileName] = optionalFileInput(fileName,'*.txt','Please select a STORM file to open:');

if filePath  == 0
    %If the user clicked cancel on file dialogue
    stormData = [];
    return
end

ip = inputParser;
ip.addParamValue('Verbose',true,@islogical);
ip.parse(varargin{:});
p = ip.Results;

%% ------ File type determination Import ---- %%

%Determine which format the specified file is
formatName = '';
for j = 1:nReader
            
    if testFun{j}([filePath fileName])
        
        formatName = readerName{j};
        iFormat = j;
        break
        
    end        
    
end

if isempty(formatName)
    error('Could not recognize file format!')
elseif p.Verbose
    disp(['File recognized as ' formatName ' STORM file.']);
end
   
%% ----- File reading and conversion ----- %%

%Read the file using the selected reader
if p.Verbose;disp('Reading STORM text file...');end
stormData = readerFun{iFormat}([filePath fileName]);

%Convert to a common data structure
if p.Verbose;disp('Converting...');end

stormData = convertStormData(stormData,formatName);

if p.Verbose;disp('Done!');end
