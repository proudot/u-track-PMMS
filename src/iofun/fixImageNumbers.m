function fixImageNumbers(directory)
%FIXIMAGENUMBERS pads zeros on metamorph-style image numbering so ls/dir return them correctly
%
% fixImageNumbers(directory)
% 
% This function goes through every numbered image file in the specified
% directory and renames files which have a metamorph-style ending like _t1,
% _t2 etc. To _t01, _t02 so that dir/ls returns them in the correct order.
% 
% Hunter Elliott, 11/2009
%

stackFiles = imDir(directory);

nStack = length(stackFiles);

nDig = floor(log10(nStack))+1;
fString = ['%0' num2str(nDig) '.0f'];


if nStack > 0
    
    disp(['Renaming ' num2str(nStack) ' image files...'])
    
    for i = 1:nStack
        
        %Get the index of the last t
        iT = max(regexp(stackFiles(i).name,'t'));
        
        %And of the file extension
        iLFS = max(regexp(stackFiles(i).name,'\.'));
        
        %Convert string to a frame number
        iFrame = str2double(stackFiles(i).name(iT+1:iLFS-1));
        
        %Make sure we have a valid number to avoid destroying data!
        if isnan(iFrame) || ~isfinite(iFrame)
            error(['Frame number not recognized! "' stackFiles(i).name(iT+1:iLFS-1) '" is not a valid frame number! Check image name format!'])
        end
        
        %Make the new file name
        newName = [stackFiles(i).name(1:iT) num2str(iFrame,fString) stackFiles(i).name(iLFS:end)];
        
        %If the file name has changed, rename it.
        if ~strcmp(newName,stackFiles(i).name)
            movefile([directory filesep stackFiles(i).name],...
                [directory filesep newName]);                
        end
    end
else 
    disp('No image files found in specified directory...')
    return
end

disp('Finished!')




