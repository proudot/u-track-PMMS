function status = isFileZeissSTORM(fileName)

%Reference header names for comparison
zeissHeaders = {'Index'    'First Frame'    'Number Frames'    'Frames Missing'    'Position X [nm]',...
                'Position Y [nm]'    'Precision [nm]'    'Number Photons'    'Background variance',...
                'Chi square'    'PSF width [nm]'    'Channel'    'Z Slice'};

status = false;
fileHeaders = {};

if exist(fileName,'file')
       
    %Attempt to read the header row from the file
    delimiter = '\t';
    startRow = 1;
    endRow = 1;
    formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
    
    try
        fileID = fopen(fileName,'r');

        textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
        fileHeaders = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);        
        fileHeaders = [fileHeaders{:}];
        
    catch                
        
    end
        
    fclose(fileID);        
    %Compare the header to the expected zeiss header
    status = isequal(zeissHeaders,fileHeaders);
        
end


       
