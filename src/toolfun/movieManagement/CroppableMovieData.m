classdef CroppableMovieData < MovieData
    %CroppableMovieData MovieData subclass that can be cropped
    
    properties
        cropDataFile_
    end
    
    methods
        function obj = CroppableMovieData(varargin)
            % cropDataFile_ can be set as a parameter
            obj@MovieData(varargin{:});
            if(isempty(obj.cropDataFile_))
                obj.cropDataFile_ = [obj.outputDirectory_ filesep 'cropReader.mat'];
            end
        end
        function r = initReader(obj)
            r = initReader@MovieData(obj);
            % if cropDataFile_ exists, then create a CropReader that
            % will crop the images from the previous reader
            if(exist(obj.cropDataFile_,'file'))
                S = load(obj.cropDataFile_);
                r = CropReader(r,S.positions{end});
            end
        end
        function setReader(obj, r)
            obj.setReader@MovieData(r);
            % On setting a reader, modify the dimensions to match the
            % dimensions of the new reader in terms of XYZT
            if(~isempty(r))
                % could be empty if trying to reinitialize the Reader
                obj.imSize_ = [ r.getSizeY() r.getSizeX() ];
                obj.nFrames_ = r.getSizeT();
                obj.zSize_ = r.getSizeZ();
            end
        end
        function crop(obj, position)
            % position is [x y w h] relative to the current image
            % dimensions
            
            if(nargin < 2)
                % if no position given, run movieViewer and use imrect to
                % acquire it
                h = movieViewer(obj);
                hr = imrect;
                wait(hr);
                position = round(getPosition(hr));
            end
            
            % find all the crop readers in the ProxyReader chain so that we
            % can reference the current position to the original movie
            cropReaders = CropReader.findCropReaders(obj.getReader());
            for cri = 1:length(cropReaders)
                r = cropReaders{cri};
                % shift by previous crop offset if it exists
                position = position + [r.position([1 2]) 0 0];
            end
            
            try
                if(exist(obj.cropDataFile_,'file'))
                    S = load(obj.cropDataFile_);
                else
                    S.positions = {};
                end

                S.positions{end+1} = position;
                save(obj.cropDataFile_,'-struct','S');
            catch err
                disp(['Could not save crop area to ' obj.cropDataFile_]);
                disp(err);
            end
            if(nargin < 2)
                delete(hr);
                close(h);
            end
            % reinitialize the reader
            obj.setReader(obj.initReader);
        end
    end
    
end

