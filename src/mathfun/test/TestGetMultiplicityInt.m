classdef TestGetMultiplicityInt < TestGetMultiplicity
    methods
        function self = TestGetMultiplicityInt(varargin)
            self = self@TestGetMultiplicity(varargin{:});
            self.fcn = @getMultiplicityInt;
        end
        function setUp(self)
            self.setUp@TestGetMultiplicity();
            self.M = uint16(self.M);
            self.udata = uint16(self.udata);
            self.sdata = uint16(self.sdata);
        end
        function testDoubleInput(self)
            [rep, udata, sdata] = self.fcn(double(self.M));
            assertEqual( class(udata) , 'double');
            assertEqual( class(sdata) , 'double');
            self.checkOutput(rep,uint16(udata),uint16(sdata));
        end
        function testUint16Input(self)
            [rep, udata, sdata] = self.fcn(uint16(self.M));
            assertEqual( class(udata) , 'uint16');
            assertEqual( class(sdata) , 'uint16');
            self.checkOutput(rep,udata,sdata);
        end
        function testNonFiniteInput(self)
            M = [NaN; Inf; double(self.M); Inf; NaN]; 
            [rep, udata, sdata] = self.fcn(M);
            self.checkOutput(rep,uint16(udata),uint16(sdata));
        end
        function testUint16MatrixInput(self)
            [rep, udata, sdata] = self.fcn(uint16(reshape(self.M,1e3,1e3)));
            assertEqual( class(udata) , 'uint16');
            assertEqual( class(sdata) , 'uint16');
            self.checkOutput(rep,uint16(udata),uint16(sdata));
        end

    end
end
