classdef TrackingPackage < Package
    % An abstract class for a geeneric Tracking Package
    
    methods
        function obj = TrackingPackage(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'TrackingPackage'];
            end
                 
            % Call the superclass constructor
            obj = obj@Package(super_args{:});        
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'U-Track';
        end 
        function m = getDependencyMatrix(i,j)   
            m = [0 0 0;  %1 DetectionProcess
                 1 0 0;  %2 TrackingProcess
                 1 1 0;];%3 PostTrackingProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = trackingPackageGUI(varargin{:});
        end
        function classes = getProcessClassNames(index)
            classes = {
                'DetectionProcess',...
                'TrackingProcess',...
                'PostTrackingProcess'};
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        function objects = getConcretePackages()
            objects(1).name = 'Single particles';
            objects(1).packageConstr = @UTrackPackage;
            objects(2).name = 'Microtubules plus-ends';
            objects(2).packageConstr = @PlusTipTrackerPackage;                        
            objects(3).name = 'Nuclei';
            objects(3).packageConstr = @NucleiTrackingPackage;
        end
        
    end 
end