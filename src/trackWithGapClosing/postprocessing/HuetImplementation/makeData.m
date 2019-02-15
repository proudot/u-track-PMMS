function [trajMatrix, trajectories] = makeData (trajectoryCharacteristics, dimension, timeStep,trajNum);

% MAKEDATA function that creates trajectories with different stages of
%          movement and converts them into the right inout format for the
%          analyzing function transBehaviorAnalysis
%
% SYNOPSIS [trajectoryMatrix, traj] = makeData (trajectoryCharacteristics,
%                                     dimension, timeStep,TrajNum)
%
% INPUT  n-by-7 matrix with columns:
%
%            startTime = start timepoint for simulated movement[time units]
%
%            endTime   = end timepoint for simulated movement[time units]
%
%            diffCoeffitient = diffusion coeffitient for specified movement
%
%            cRadius   = Confinement radius [space units]. Optional. 
%                        Needed only in case of constrained motion.
% 
%            driftVel : Drift velocity in spherical polar coordinates. 
%                       1D: v; 2D: v,theta; 3D: v,theta,phi.
%               driftV      = speed in [(space units)/(unit time)]
%
%               driftTheta  = in degrees (only with 2D)
% 
%               driftPhi    = in degrees (only with 3D)
%
%        dimension = Dimension of space: 1, 2 or 3.
%
%        timeStep  = Simulation time step [time units].
%                    (Optional. Default: 0.1/diffConst.)
%
%
%Output trackMatrix: n-by-4 matrix with columms:
%
%       timePoints = time points of simulated movement[time units]
%
%       x          = x coordinates
%
%       y          = y coordinates
%
%       z          = z coordinates
%
% created by gp 03/29/07


for i= 1:trajNum
    trajectoryMatrix = [];
    
[traj] = brownianMotionGunnarEdit(trajectoryCharacteristics,dimension, timeStep);
% figure;
% hold on;

%plot3(traj(:,2),traj(:,3),traj(:,4));
plot(traj(:,2),traj(:,3),'g');
minDiffX = min(abs(diff(traj(:,2))));

minDiffY = min(abs(diff(traj(:,3))));

dX = minDiffX / 2;
dY = minDiffY / 2;

trajectoryMatrix(:,1) = traj(:,2);
trajectoryMatrix(:,2) = traj(:,3);
trajectoryMatrix(:,3:4) = 0;
trajectoryMatrix(:,5) = minDiffX;
trajectoryMatrix(:,6) = minDiffY;
trajectoryMatrix(:,7:8) = 0;
[rowSizeTraj, columnSizeTraj] = size(traj);
trajectoryMatrix = trajectoryMatrix';
trajMatrix(i,:) = reshape(trajectoryMatrix, 1, rowSizeTraj*8);
trajectories(:,:,i)=traj;
end
end
