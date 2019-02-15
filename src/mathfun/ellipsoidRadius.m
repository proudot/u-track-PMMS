function radius = ellipsoidRadius(vectors,majorAxis)
%ellipsoidRadius calculates the radius of an n-dimensional ellipsoid in the direction of the vectors
%
%SYNOPSIS radius = ellipsoidRadius(vectors,majorAxis)
%
%INPUT      vectors: m-by-n matrix, m vectors with n dimensions indicating the
%                       directions in which the radius has to be measured
%           majorAxis: m-by-n or 1-by-n matrix containing the major axis of the
%                       ellipsoid(s) (e.g. for a circle with radius 1: [1,1])
%
%OUTPUT     radius: m-by-1 matrix with lengths of the radius in the specified
%                       directions
%
%c: 03/05/27 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%test input
vSize = size(vectors);
mASize = size(majorAxis);

if length(vSize)>2|length(mASize)>2
    error('input has to be 2D-arrays')
end

if vSize(2) ~= mASize(2)
    error('both inputs have to be of the same dimension (# of cols)')
end

if ~( vSize(1)==mASize(1) | mASize(1)==1 )
    error('specify either one ellipsoid for all inputs or one ellipsoid for every input vector')
end

if any(majorAxis==0)
    warning('at least one major axis has length 0!')
end

%norm vectors
[dummy, unitVectors] = normList(vectors);

%make majorAxis the right size
if mASize(1)==1 & vSize(1)~=1
    majorAxis = ones(vSize(1),1)*majorAxis;
end

%calculate radius: radius = 1/sqrt(vx^2/rx^2 + vy^2/ry^2 + vz^2/rz^2); 
%v = unit vector, r = major axis
sumOfSquares = sum( (unitVectors.^2) ./ (majorAxis.^2),2 );
radius = sqrt(sumOfSquares).^-1;