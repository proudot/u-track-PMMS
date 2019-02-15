function rvec=rot3vec(vec,direction,alpha,origin)
%ROT3VEC rotate point set about direction vector through origin and angle alpha
%
% SYNOPSIS rvec=rot3vec(vec,direction,alpha,origin)
%
% INPUT vec         : 3D point set [x(:) y(:) z(:)]
%       direction   : 3D vector
%       alpha       : angle in degree
%       origin      : 3D vector center of rotation 
%
% OUTPUT rvec       : rotated 3D point set [rx(:) ry(:) rz(:)]

% c: 21/09/01	dT


if nargin <4
    origin=[0 0 0];
end;
R = rot3(direction,alpha);

newxyz = [vec(:,1)-origin(1), vec(:,2)-origin(2), vec(:,3)-origin(3)];
newxyz = newxyz*R;
rvec = [newxyz(:,1)+origin(1),newxyz(:,2)+origin(2), newxyz(:,3)+origin(3)];



