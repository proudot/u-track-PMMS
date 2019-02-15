function status = isFileNikonSTORM(fileName)

%Reference header names for comparison
nikonHeaders = {'Channel Name'    'X'    'Y'    'Xc'    'Yc'    'Height'    'Area'    'Width'    'Phi'    'Ax',...
                'BG'    'I'    'Frame'    'Length'    'Link'    'Valid'    'Z'    'Zc'};

status = false;
fileHeaders = {};

if exist(fileName,'file')
       
    %Attempt to read the header row from the file
    delimiter = '\t';
    startRow = 1;
    endRow = 1;
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
    
    try
        fileID = fopen(fileName,'r');
        fileHeaders = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
        fileHeaders = [fileHeaders{:}];
          
    catch                
        
    end
        
    fclose(fileID);        
    %Compare the header to the expected nikon header
    status = isequal(nikonHeaders,fileHeaders);
        
end