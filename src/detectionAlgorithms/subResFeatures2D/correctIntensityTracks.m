function trackFinal = correctIntensityTracks(tracksFinal,correctionImage)
%CORRECTINTENSITYTRACKS corrects object intensities in tracksFinal due to uneven illumination
%
%SYNOPSIS trackFinal = correctIntensityTracks(tracksFinal,correctionImage)
%
%INPUT  
%       tracksFinal  : Tracks information, e.g. output of
%                      trackCloseGapsKalmanSparse.
%       correctionImage: Image to use for correction.
%
%OUTPUT tracksFinal  : Same as input, but with intensities replaced with
%                      corrected ones.
%
%Khuloud Jaqaman, August 2014

%% Input

%find number of frames in movie
numFrames = length(movieInfo);

%find image size
correctionImage = double(correctionImage);
imSize = size(correctionImage);

%check if raw intensities are input
if isfield(movieInfo,'rawInt')
    fieldsNew = {'amp','rawInt','rawIntMinusBg'};
    fieldsOld = {'ampUncorr','rawIntUncorr','rawIntMinusBgUncorr'};
else
    fieldsNew = {'amp'};
    fieldsOld = {'ampUncorr'};
end
numField = length(fieldsNew);

%% Correction

%normalize correction image so that average is 1
meanInt = mean(correctionImage(:));
correctionImage = correctionImage / meanInt;

%go over all frames
for iFrame = 1 : numFrames
   
    %move uncorrected amplitudes/intensities to new field
    for iField = 1 : numField
        movieInfo(iFrame).(fieldsOld{iField}) = movieInfo(iFrame).(fieldsNew{iField});
    end
    
    if ~isempty(movieInfo(iFrame).xCoord)
        
        %object information
        xCoord = round(movieInfo(iFrame).xCoord(:,1));
        yCoord = round(movieInfo(iFrame).yCoord(:,1));
        for iField = 1 : numField
            amp(:,iField) = movieInfo(iFrame).(fieldsNew{iField})(:,1);
        end
        
        %linear index
        linIdx = sub2ind(imSize,yCoord,xCoord);
        
        %correction image value at linIdx
        correctionVal = correctionImage(linIdx);
        
        %correct amplitude/intensity
        correctionVal = repmat(correctionVal,1,numField);
        amp = amp ./ correctionVal;
        for iField = 1 : numField
            movieInfo(iFrame).(fieldsNew{iField})(:,1) = amp(:,iField);
        end
        
    end
    
end

%%%%% ~~ the end ~~ %%%%%
