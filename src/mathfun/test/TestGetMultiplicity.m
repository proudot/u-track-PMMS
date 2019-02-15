classdef TestGetMultiplicity < TestCase
    properties
        M
        fcn
        rep
        udata
        sdata
    end
    methods
        function self = TestGetMultiplicity(varargin)
            self = self@TestCase(varargin{:});
            self.fcn = @getMultiplicity;
        end
        function setUp(self)
            self.M = randi(intmax('uint16'),1e6,1);
            self.udata = unique(self.M)';
            self.sdata = sort(self.M)';
            self.rep = histc(self.M, self.udata)';
        end
        function checkOutput(self,rep,udata,sdata)
            assertEqual(rep, self.rep);
            if(nargin > 2)
                assertEqual(udata, self.udata);
            end
            if(nargin > 3)
                assertEqual(sdata, self.sdata);
            end
        end
	function testEmpty(self)
		[rep, udata, sdata] = self.fcn([]); 
		assert(isempty(rep));
		assert(isempty(udata));
		assert(isempty(sdata));
		rep = self.fcn([]);
		assert(isempty(rep));
	end
        function testSingleOutput(self)
            self.checkOutput(self.fcn(self.M));
        end
        function testDoubleOutput(self)
            [rep, udata] = self.fcn(self.M);
            self.checkOutput(rep,udata);
        end
        function testTripleOutput(self)
            [rep, udata, sdata] = self.fcn(self.M);
            self.checkOutput(rep,udata,sdata);
        end
        function testMatrixInput(self)
            [rep, udata, sdata] = self.fcn(reshape(self.M,1e3,1e3));
            self.checkOutput(rep,udata,sdata);
        end
        function testNonFiniteInput(self)
            M = [NaN; Inf; self.M; Inf; NaN]; 
            [rep, udata, sdata] = self.fcn(M);
            self.checkOutput(rep,udata,sdata);
        end
        function testUint16Input(self)
            [rep, udata, sdata] = self.fcn(uint16(self.M));
            assertEqual( class(udata) , 'uint16');
            assertEqual( class(sdata) , 'uint16');
            self.checkOutput(rep,double(udata),double(sdata));
        end
        function testUint16ExtremesInput(self)
            % type uint16
            M = [0 1 intmax('uint16')];
            [rep, udata, sdata] = self.fcn(M);
            assertEqual( rep, [ 1 1 1]);
            assertEqual( udata, M);
            assertEqual( sdata, M);
        end
        function testUint16MatrixInput(self)
            [rep, udata, sdata] = self.fcn(uint16(reshape(self.M,1e3,1e3)));
            assertEqual( class(udata) , 'uint16');
            assertEqual( class(sdata) , 'uint16');
            self.checkOutput(rep,double(udata),double(sdata));
        end
        function tearDown(self)
            self.M = [];
        end
    end
end
