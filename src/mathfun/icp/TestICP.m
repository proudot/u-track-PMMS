classdef TestICP <TestCase
    % This function test the Iterative Closest Point algorithm (ICP).
    
    % Adapated from testComputeICP.m
    
    properties
        tol = 1e-9; % ending condition of ICP
        prec = 10; % number of significant digit
        maxIter = 10;
    end
    methods
        function self = TestICP(name)
            self = self@TestCase(name);
        end
        
        function testIdentity(self)
            X1 = rand(100, 3);
            
            [T R] = computeICP(X1, X1, self.maxIter, self.tol);
            T = round(T * 10^self.prec) / 10^self.prec;
            R = round(R * 10^self.prec) / 10^self.prec;
            assertTrue(all(T == 0) && all(all(R == eye(3))))
        end
        
        
        function testTranslation(self)
            X1 = rand(100, 3);
            X2 = X1;
            shift = [1e-9 0 0];
            X1 = bsxfun(@plus, X1, shift);
            
            [T R] = computeICP(X1, X2, self.maxIter, self.tol);
            T = round(T * 10^self.prec) / 10^self.prec;
            R = round(R * 10^self.prec) / 10^self.prec;
            assertTrue(all(T == shift') && all(all(R == eye(3))));
        end
        
        function testRotation(self)
            X1 = rand(100, 3);
            theta = (.5 - rand(1)) * pi;
            ct = cos(theta);
            st = sin(theta);
            Rref = [ct st 0; -st ct 0; 0 0 1];
            X2 = X1;
            c1 = mean(X1);
            X1 = bsxfun(@minus, X1, c1);
            for i = 1:size(X1, 1)
                X1(i, :) = X1(i, :) * Rref;
            end
            X1 = bsxfun(@plus, X1, c1);
            
            [T R] = computeICP(X1, X2, self.maxIter, self.tol);
            T = round(T * 10^self.prec) / 10^self.prec;
            R = round(R * 10^self.prec) / 10^self.prec;
            assertTrue(all(T == 0) && abs(acos(R(1))^2 - theta^2) < 1e-9 && R(end) == 1)
        end
        
        function testTranslationRotaion(self)
            X1 = rand(100, 3);
            shift = [1e-9 0 0];
            theta = (.5 - rand(1)) * pi;
            ct = cos(theta);
            st = sin(theta);
            Rref = [ct st 0; -st ct 0; 0 0 1];
            X2 = X1;
            c1 = mean(X1);
            X1 = bsxfun(@minus, X1, c1);
            for i = 1:size(X1, 1)
                X1(i, :) = X1(i, :) * Rref;
            end
            X1 = bsxfun(@plus, X1, c1 + shift);
            
            [T R] = computeICP(X1, X2, self.maxIter, self.tol);
            T = round(T * 10^self.prec) / 10^self.prec;
            R = round(R * 10^self.prec) / 10^self.prec;
            assertTrue(all(T == shift') && abs(acos(R(1))^2 - theta^2) < 1e-9 && R(end) == 1)
        end
    end
end

