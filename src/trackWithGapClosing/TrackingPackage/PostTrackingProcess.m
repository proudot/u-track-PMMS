classdef PostTrackingProcess < DataProcessingProcess
    % An abstract class for all post-tracking processes 
    
    % Sebastien Besson Jul 2012
    
    methods(Access = public)
        
        function obj = PostTrackingProcess(owner, name, funName, funParams )
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;
            end
            if nargin > 3
                obj.funParams_ = funParams;
            end
        end
        
    end
    methods(Static)
        function name = getName()
            name = 'Track analysis';
        end
        function h = GUI()
            h = @abstractProcessGUI;
        end
        function procClasses = getConcreteClasses()
            procClasses = ...
                {'MotionAnalysisProcess';
                'CometPostTrackingProcess';                
                'TrackGroupingProcess'};
        end
    end
end