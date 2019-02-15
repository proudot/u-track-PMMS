classdef ExternalProcess < Process
    % A concrete class
    
    methods(Access = public)
        
        function obj = ExternalProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('name',ExternalProcess.getName(),@ischar);
            ip.parse(owner,varargin{:});
            
            % Constructor of the DummyDetectionProcess
            super_args{1} = owner;
            super_args{2} = ip.Results.name;
            obj = obj@Process(super_args{:});
            obj.funName_ = @(x) x;
            obj.funParams_ = ExternalProcess.getDefaultParams(varargin);
            
        end
    end
    methods (Static)
        function name = getName()
            name = 'Dummy process';
        end
        function funParams = getDefaultParams(varargin)
            funParams = struct();
        end
        
    end
end
