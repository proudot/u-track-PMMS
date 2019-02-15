function R=rot3(direction,alpha,degSwitch)
%ROT3 create a rotation matrix from direction vector and angle
%
% SYNOPSIS R=rot3(direction,alpha,degSwitch)
%
% INPUT direction   : 3D vector
%       alpha       : angle in degree or radian (see degSwitch!)
%       degSwitch   : 'deg' (default) or 'rad', depending on alpha input
%
% OUTPUT
%       R           : rotation matrix (clockwise rotation; inverse rotAxis-direction for counterclockwise rotation)

% c: 21/09/01	dT

%check input
if nargin==2|isempty(degSwitch)
    degSwitch='deg';
elseif ~(strcmp(degSwitch,'deg')|strcmp(degSwitch,'rad'))
    %disp warning, switch to rad
    disp(['warning: rot3 called with unknown degSwitch ''',degSwitch,'''; switching to ''rad''']);
    degSwitch='rad';
end

% normalize  direction
if all(direction ==0)
    R=eye(3);
    return;
end;
direction=direction(:)/norm(direction);
% switch to radians if necessary
switch degSwitch
    case 'deg'
        alph = alpha*pi/180;
    case 'rad'
        alph=alpha;
end

cosa = cos(alph);
sina = sin(alph);
vera = 1 - cosa;
x = direction(1);
y = direction(2);
z = direction(3);
R = [cosa+x^2*vera x*y*vera-z*sina x*z*vera+y*sina; ...
     x*y*vera+z*sina cosa+y^2*vera y*z*vera-x*sina; ...
     x*z*vera-y*sina y*z*vera+x*sina cosa+z^2*vera]';
