
% Lower and upper bound test of the computed Bezier curve length 

clear all; clc;

nTests = 100000;
scalingFactor = 10;
curveDeg = 3;

for i=1:nTests
    
    cP = rand(curveDeg+1,3)*scalingFactor;
    
    lowLenBound = norm(cP(1,:)-cP(end,:));
    upLenBound = sum(sqrt(sum((cP(1:end-1,:)-cP(2:end,:)).^2,2)));
    
    len = lengthBezier(cP);
    
    % fprintf('%f %f %f\n',lowLenBound,len,upLenBound);
    
    if len < lowLenBound
        disp('Computed curve is too short!');
    elseif len > upLenBound
        disp('Computed curve is too long!');
    else
        % disp('OK!');
    end
    
end