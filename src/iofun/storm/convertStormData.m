function stormData = convertStormData(stormData,formatName)
%CONVERTSTORMDATA converts storm data from different proprietary formats to a common structure
%
% stormData = convertStormData(stormData,formatName)
%
% UNDER CONSTRUCTION!! Currently only standardizes the X and Y
% coordinates....
%
%   Input:
%       A STORM data structure as produced by the dedicated reader function
%       for that format (see e.g. readNikonSTORM, readZeissSTORM)
%
%   Output:
%       stormData - a structure with standardized field names, AND the
%       non-standard fields from each format
%

%Hunter Elliott
%1/2015

switch formatName
    
    
    case 'Nikon'
        
        %Just waste memory who gives a fuck. Store the un-corrected X and Y
        %before replacing
        tmpX = stormData.X;
        tmpY = stormData.Y;
        
        %Replace them with the drift-corrected coord
        stormData.X = stormData.Xc;
        stormData = rmfield(stormData,'Xc');
        stormData.X_uncorrected = tmpX;%Store the old coord specifying they are uncorrected
        
        stormData.Y = stormData.Yc;
        stormData = rmfield(stormData,'Yc');
        stormData.Y_uncorrected = tmpY;%Store the old coord specifying they are uncorrected
        
        %TEMP - get frame numbers, other fields and convert!!
        
        
    case 'Zeiss'
        
        %"rename" the positions by copy & delete.
        stormData.X = stormData.PositionXnm;
        stormData = rmfield(stormData,'PositionXnm');
        stormData.Y = stormData.PositionYnm;
        stormData = rmfield(stormData,'PositionYnm');
        
        %TEMP - get frame numbers etc and convert!
        
    otherwise
        
        error('unrecognized format name!')
end
        
        
        



