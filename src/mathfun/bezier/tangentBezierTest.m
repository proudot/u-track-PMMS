
clear all; clc;

% Initialize a curve and a few nodes
P = [[1,2];[1,5];[4,2];[7,8]];
t = [0;0.2;0.4;0.9;1];

% Compute tangents
[T,normalT] = tangentBezier(P,t);
T = T./repmat(sqrt(sum(T.^2,2)),1,2);
Y = renderBezier(P,t);

% Render curve
ts = linspace(0,1,50)';
Ys = renderBezier(P,ts);
figure(1);
plot(Ys(:,1),Ys(:,2));

% Draw Tangents
hold on;
for i=1:numel(t)
    if T(i) ~= normalT(i)
        disp('Oops!');
    end
    plot([Y(i,1)-T(i,1),Y(i,1)+T(i,1)],[Y(i,2)-T(i,2),Y(i,2)+T(i,2)],'r-')
end
hold off;