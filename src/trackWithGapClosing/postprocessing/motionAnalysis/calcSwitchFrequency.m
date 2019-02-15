function switchFreq = calcSwitchFrequency(trajectory)
%CALCSWITCHFREQUENCY calculate the directional switching frequency in a trajectory
%
%SYNPOSIS
%
%INPUT
%
%OUTPUT
%
%Khulud Jaqaman, March 2010

%% Output

switchFreq = [];

%% Input

if nargin < 1
    disp('Please enter trajectory')
    return
end

%% Calculation

%calculate displacements
dispTraj = diff(trajectory);

%calculate the dot product between consecutive displacements
dispDotProduct = dot(dispTraj(1:end-1,:),dispTraj(2:end,:),2);

%calculate switching frequency
switchFreq = length(find(dispDotProduct<0)) / length(dispDotProduct);

%% ~~~ the end ~~~