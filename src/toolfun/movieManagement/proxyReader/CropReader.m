classdef CropReader < ProxyReader
    %CropReader Crops an image in X and Y using imcrop
    
    properties
        position;
    end
    
    methods
        function obj = CropReader(varargin)
            if(nargin < 1)
                % no parameters, return nothing to allow for extension
                return;
            end
            if(nargin > 1)
                obj.position = round(varargin{2});
            end
            % at least one parameter
            obj.setReader(varargin{1});
        end
        function oldReader = setReader(obj,reader)
            oldReader = obj.setReader@ProxyReader(reader);
            if(isempty(obj.position))
                hfig = figure;
                imshow(reader.loadImage(1,1),[]);
                h = imrect;
                wait(h);
                obj.position = getPosition(h);
                delete(h);
                close(hfig);
            end
        end
        function sizeY =  getSizeY(obj)
            sizeY = obj.position(4)+1;
        end
        function sizeX = getSizeX(obj)
            sizeX = obj.position(3)+1;
        end
        function I = loadImage(obj,varargin)
            I = obj.loadImage@ProxyReader(varargin{:});
            I = imcrop(I,obj.position);
        end
        function I = loadStack(obj, c, t, varargin)
            I = loadStack@Reader(obj, c, t, varargin{:});
        end
    end
    methods ( Access = protected )
        function I = loadImage_(obj,varargin)
            I = obj.loadImage_@ProxyReader(varargin{:});
            I = imcrop(I,obj.position);
        end
        function I = loadStack_(obj, c, t, z)
            % do not proxy, use default implementation
            I = loadStack_@Reader(obj, c, t, z);
        end
    end
    methods ( Static )
        function cropReaders = findCropReaders(reader)
            cropReaders = [];
            if(isa(reader,'ProxyReader'))
                readers = reader.findProxies();
                isCropReader = cellfun(@(r) isa(r,'CropReader'),readers);
                cropReaders = readers(isCropReader);
            end
        end
        function [reader] = apply(movieData,position)
            % Applies the CropReader to movieData
            % INPUT
            % movieData - MovieData object instance
            % position - (optional) [x y w h] position
            % should be replaced (default: true)
            %
            % OUTPUT
            % reader - a CropReader instance
            if(nargin < 2)
                h = movieViewer(movieData);
                hr = imrect;
                wait(hr);
                position = round(getPosition(hr));
            end
            reader = CropReader(movieData.reader,position);
            movieData.setReader(reader);
            file = [ movieData.outputDirectory_ filesep 'cropReader.mat'];
            try
                if(exist(file,'file'))
                    matData = load(file);
                else
                    matData.positions = {};
                end
                matData.positions{end+1} = position;
                save(file,'-struct','matData');
            catch err
                delete(hr);
                close(h);
                disp(['Could not save crop area to ' file]);
                disp(err);
                return;
            end
            delete(hr);
            close(h);
        end
    end
end

