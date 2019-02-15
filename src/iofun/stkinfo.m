function info = stkinfo(fileName)
%STKINFO reads and returns header information from metamorph STK TIFF file
% 
% info = stkinfo(fileName)
% 
% This function reads various imaging parameters from the header of the
% specified metamorph .STK file and returns them. 
%
% Redundancy with other functions: 
% Although the tiffread.m function (and others) can also return information
% from the STK header, the difference is that this function does so without
% reading the actual image data from the stack, making it faster. It also
% returns some additional information which is not returned by tiffread.m
% or imfinfo.m.
% 
% This function borrows very heavily from tiffread.m so thanks a bunch to
% Francois Nedelec et al.
%
% *NOTE* This function assumes that the image information is the same for
% every image in the stack - the information is only read from the first
% image in the stack.
% 
% Input:
% 
%   fileName - A character array specifiying the name of the .STK file to
%   read information from. This input should include the file extension.
%   Optional. If not input, the user will be asked to select a file.
% 
% Output:
% 
%   info - A structure containing all the image information. The field
%   names should be self-explanatory, and any additional user-entered
%   information is stored in the Notes field.
% 
% 
% Hunter Elliott
% 4/2011

%% ---------------- Input ------------------ %%

if nargin < 1 || isempty(fileName)    
    [fName,path] = uigetfile({'*.STK','*.stk'},'Select an image:');
    if fName == 0
        info = [];
        return
    else
        fileName = [path filesep fName];
    end
elseif ~exist(fileName,'file')
    error(['The file ' fileName ' does not exist!'])
end


%% ----------------- Init ------------------- %%


%Open the file for reading
[TIF.file,sysMess] = fopen(fileName,'r','l');
if TIF.file == -1
    error(['The file ' fileName ' could not be opened! : ' sysMess]);
end

%Determine bite ordering and make sure this is a valid STK/TIFF file. Store
%the Byte Ordering String for later use.
byteOrder = fread(TIF.file,2,'*char');
if strcmp(byteOrder', 'II')
    TIF.BOS = 'ieee-le';                                % Intel little-endian format
elseif ( strcmp(byteOrder','MM') )
    TIF.BOS = 'ieee-be';
else
    error(['The file ' fileName ' is not a valid STK file!']);
end

%The number 42 ("the answer") further identifies the file as a .tif file.
%Make sure that this is present
tiffId = fread(TIF.file,1,'uint16', TIF.BOS);
if (tiffId ~= 42)
    error(['The file ' fileName ' is not a valid STK file!']);
end

%Now we find the location of the first Image File Directory (IFD) in the file
TIF.imgPos = fread(TIF.file, 1, 'uint32', TIF.BOS);

status = fseek(TIF.file,TIF.imgPos,-1);
if status == -1
    error('Error reading file: invalid file offset (error on fseek)!')
end

%Check the number of entries/fields in this IFD
nEntries = fread(TIF.file,1,'uint16',TIF.BOS);

%Go through each entry, and read and store all supported types. 
for i = 1:nEntries

    %Get the current position in the file
    filePos = ftell(TIF.file);

    %Read the tag for this entry.The tag numbers for each entry are
    %defined by the adobe TIF standard and additional tags are
    %defined by metamorph.
    TIF.entryTag = fread(TIF.file,1,'uint16',TIF.BOS);

    %Read the actual entry value
    entry = readIFDentry(TIF);

    %Store this in an appropriately named field, based on the entry tag
    %number. This switch statement is largely copied from tiffread.m,
    %but modifed: More fields are returned, and the field names use the
    %TIFF-standard 6.0 defined names rather than the tiffread names.
    switch TIF.entryTag
        case 254
            info.NewSubfiletype = entry.val;
        case 256         % image width - number of column
            info.Width          = entry.val;
        case 257         % image height - number of row
            info.Height         = entry.val;
            info.ImageLength    = entry.val;
        case 258         % BitsPerSample per sample
            info.BitsPerSample  = entry.val;
            info.BytesPerSample = info.BitsPerSample / 8;
            info.Bits           = info.BitsPerSample(1);                
        case 259         % compression
            info.CompressionType = entry.val;                
        case 262         % photometric interpretation
            info.PhotometricInterpretation = entry.val;                
        case 269
            info.DocumentName  = entry.val;
        case 270         % comments:
            info.ImageDescription = entry.val;
        case 271
            info.Make           = entry.val;
        case 273         % strip offset
            info.StripOffsets   = entry.val;
            info.StripNumber    = entry.cnt;                
        case 277         % sample_per pixel
            info.SamplesPerPixel  = entry.val;                
        case 278         % rows per strip
            info.RowsPerStrip   = entry.val;
        case 279         % strip byte counts - number of bytes in each strip after any compressio
            info.StripByteCounts= entry.val;
        case 282         % X resolution
            info.XResolution   = entry.val;
        case 283         % Y resolution
            info.YResolution   = entry.val;
        case 284         %planar configuration describe the order of RGB
            info.PlanarConfiguration = entry.val;
        case 296         % resolution unit
            info.ResolutionUnit= entry.val;
        case 305         % software
            info.Software       = entry.val;
        case 306         % datetime
            info.DateTime       = entry.val;
        case 315
            info.Artist         = entry.val;
        case 317        %predictor for compression
            info.Predictor = entry.val;
        case 320         % color map
            info.ColorMap           = entry.val;
            info.Colors         = entry.cnt/3;
        case 339
            info.SampleFormat   = entry.val;
        case 33628       %metamorph specific data for older metamorph versions
            info.OldMetamorphData   = entry.val;
        case 33629       %this tag identifies the image as a Metamorph stack!
            info.StackInfo      = entry.val;
            info.NumZPlanes = entry.cnt;
        case 33630       %metamorph stack data: wavelength
            info.Wavelength  = entry.val;
        case 33631       %metamorph stack data. Contains various other information on each slice - at some point I should try to read all this data also.
            info.MetamorphData    = entry.val;           
        otherwise

    end
    % move to next IFD entry in the file
    status = fseek(TIF.file, filePos+12, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end
end    
            
%Parse the metamorph-specific info in the ImageDescription entry
if isfield(info, 'StackInfo') && isfield(info, 'ImageDescription') ...
        && ~isempty(deblank(info.ImageDescription))
    info = parseMetamorphInfo(info);    
else
    warning('STKINFO:noMetamorph','No metamorph-specific file information found!');
end


%Close the file
fclose(TIF.file);                       
        
        

function entry = readIFDentry(TIF)
%Sub-function for reading IFD entry. Modified from tiffread.m to avoid
%global variable use.

entry.tiffType = fread(TIF.file, 1, 'uint16', TIF.BOS);
entry.cnt      = fread(TIF.file, 1, 'uint32', TIF.BOS);

[ entry.nbBytes, entry.matlabType ] = convertType(entry.tiffType);

if entry.nbBytes * entry.cnt > 4
    %next field contains an offset:
    offset = fread(TIF.file, 1, 'uint32', TIF.BOS);    
    status = fseek(TIF.file, offset, -1);
    if status == -1
        error('invalid file offset (error on fseek)');
    end

end


if TIF.entryTag == 33629 % metamorph 'rationals' for stack info
    entry.val = fread(TIF.file, 6*entry.cnt, entry.matlabType, TIF.BOS);
elseif TIF.entryTag == 33630   
    entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);    
else
    if entry.tiffType == 5
        entry.val = fread(TIF.file, 2*entry.cnt, entry.matlabType, TIF.BOS);
    else
        entry.val = fread(TIF.file, entry.cnt, entry.matlabType, TIF.BOS);
    end
end

if ( entry.tiffType == 2 );
    entry.val = char(entry.val');
end


function [nbBytes, matlabType] = convertType(tiffType)

%Sub function for converting the tiff data-type identifier numbers to
%matlab recognized class name strings.
%Copied from tiffread.m
switch (tiffType)
    case 1
        nbBytes=1;
        matlabType='uint8';
    case 2
        nbBytes=1;
        matlabType='uchar';
    case 3
        nbBytes=2;
        matlabType='uint16';
    case 4
        nbBytes=4;
        matlabType='uint32';
    case 5
        nbBytes=8;
        matlabType='uint32';
    case 7
        nbBytes=1;
        matlabType='uchar';
    case 11
        nbBytes=4;
        matlabType='float32';
    case 12
        nbBytes=8;
        matlabType='float64';
    otherwise
        error('tiff type %i not supported', tiffType)
end

function info = parseMetamorphInfo(info)
%Sub-routine for parsing the info that metamorph stores in the
%ImageDescription IFD.Modified from the routine used in tiffread.m to
%better handle notes and to combine the metamorph info with the regular
%tiff info, and to convert other metamorph-specific fields
%to human-readable format. 

inf   = regexprep(info.ImageDescription, '\r\n|\o0', '\n');
parse  = textscan(info.ImageDescription, '%s %s', 'Delimiter', ':');
tokens = parse{1};
values = parse{2};

% If the stk file has been created in imageJ, the delimiter used is '='.
if isempty(char(values))
    parse = textscan(inf, '%s%s', 'Delimiter', '=');
    tokens = parse{1};
    values = parse{2};
    
    if isempty(char(values))
        error('Invalid delimiter in metamorph stack info.');
    end
end

first = char(tokens(1,1));

k = 0;
iNote = 1;
for i=1:size(tokens,1)
    tok = char(tokens(i,1));
    val = char(values(i,1));
    %fprintf( '"%s" : "%s"\n', tok, val);
    if strcmp(tok, first)
        k = k + 1;
    end
    if strcmp(tok, 'Exposure')
        [v, ~, ~, pos] = sscanf(val, '%i');
        unit = val(pos:length(val));
        %return the exposure in milli-seconds
        switch( unit )
            case 'ms'
                info(k).Exposure = v;
            case 's'
                info(k).Exposure = v * 1000;
            otherwise
                warning('tiffread2:Unit', ['Exposure unit "',unit,'" not recognized']);
                info(k).Exposure = v;
        end
    else
        switch tok
            case 'Binning'
                % Binning: 1 x 1 -> [1 1]
                info(k).Binning = sscanf(val, '%d x %d')';
            case 'Region'
                info(k).Region = sscanf(val, '%d x %d, offset at (%d, %d)')';
            otherwise                
                if isempty(val)
                    %If there is no value, we store the full token string in a
                    %field called Note.
                    info(k).Notes{iNote} = tok;
                    iNote = iNote+1;                    
                else
                    field = genvarname(regexprep(tok, ' ', ''));
                    if strcmp(val, 'Off')
                        info(k).(field)=0;
                    elseif strcmp(val, 'On')
                        info(k).(field)=1;
                    elseif isstrprop(val,'digit')
                        info(k).(field)=str2num(val); %#ok<ST2NM>
                    else
                        info(k).(field)=val;
                    end
                end
        end
    end
end

%Now, convert the stack info into something useful
info.StackInfo = reshape(info.StackInfo,[6 info.NumZPlanes]);
%The z-spacing is stored as a rational. We assume all planes are equally
%spaced.
info.ZSpacing = double(info.StackInfo(1,1)) / double(info.StackInfo(2,1));

%Convert the wavelength to readable numbers. We assume all planes are the
%same wavelength
info.Wavelength = info.Wavelength(1) / info.Wavelength(2);


