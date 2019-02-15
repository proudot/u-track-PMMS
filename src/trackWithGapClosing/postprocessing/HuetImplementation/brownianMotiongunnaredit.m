function [traj] = brownianMotionGunnarEdit(trajectoryCharacteristics, dimension, timeStep)
%BROWNIANMOTION simulates constrained/unconstrained Brownian motion in 1, 
%               2 and 3D
%
%SYNOPSIS [trackMatrix] = brownianMotion(trajectoryCharacteristics, 
%                         dimensions, timeStep)
%
%INPUT  n-by-7 matrix with columns:
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
% created by gp 03/19/07

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traj = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determines the row and columnsize of the input matrix 
% trajectoryCharacteristics
[rowSizeInput, ColumnSizeInput] = size(trajectoryCharacteristics);

% check if drift velocity has the right dimensions
if ColumnSizeInput == 7 && dimension ~= 3 || ColumnSizeInput == 6 && dimension ~= 2 || ColumnSizeInput == 5 && dimension ~= 1;
   error('check dimensions');
end

%check if input values are positive or NaN
if any(trajectoryCharacteristics(:) <= 0) || any(sum(isnan(trajectoryCharacteristics(:,1:3))));
    error(['all input values in column 1:3 have to be positive all others',...
           'have to be positive or NaN']);
end

%check if correct number of arguments were used when function was called
if nargin < 2
    error('Incorrect number of input arguments!');
end



%check dimensionality of space
if ~any(dimension == [1 2 3])
    error('Variable "dimension" should be 1, 2 or 3!');
end

%initialize
countNumIterations = 1;

lastPosition(1,2:dimension+1) = 0;

%defines the total end time point
endEndTime = trajectoryCharacteristics(rowSizeInput, 2)/timeStep;

%initialize output vector
traj = zeros(endEndTime ,dimension+1);
% setting time points in the output vector
traj(:,1) = 1:endEndTime;


% loops through all rows of the input matrix trajectoryCharacteristics
for n = 1:rowSizeInput



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % assignment of variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        startTime = trajectoryCharacteristics(n,1);
        endTime = trajectoryCharacteristics(n,2);
        diffCoeffitient = trajectoryCharacteristics(n,3);
        cRadius = trajectoryCharacteristics(n,4);
        driftV = trajectoryCharacteristics(n,5);
        if dimension > 1
            driftTheta =  trajectoryCharacteristics(n,6);
        end
        if dimension > 2
            driftPhi = trajectoryCharacteristics(n,7);
        end

        %check simulation time step
        if nargin < 3 || isempty(timeStep) %if "timeStep" is not provided by user
           timeStep = 0.1/diffCoeffitient; %use default
        elseif timeStep <= 0
           error('"timeStep" should be positive!');
        end

        %check if motion is constrained
        if isnan(cRadius)
            constrained = 0;
        else
            constrained = 1;
        end

        %check if motion is directed
        if isnan(driftV)
            directed = 0;
        else
            directed =1;
        end

        %     %check drift velocity
        %     if nargin < 7 || isempty(driftVel) %if "driftVel" is not provided by user
        %        driftVelCart = zeros(1,dimension); %use default
        %     else

        %write drift velocity in Cartesian coordinates
        switch dimension
            case 1
                driftVelCart = driftV;
            case 2
                theta = driftTheta*pi/180;
                driftVelCart = driftV*[cos(theta) sin(theta)];
            case 3
                theta = driftTheta*pi/180;
                phi = driftPhi*pi/180;
                driftVelCart = driftV*[cos(theta)*sin(phi) ...
                    sin(theta)*sin(phi) cos(phi)];
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Trajectory generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points
%                 
%                 if i ==1
%                     i = 2;
%                 end
        %determine number of iterations to perform
        numIterations = ceil((endTime-startTime+1)/timeStep);

        %assuming that the step taken in a timeStep is normally distributed with mean
        %zero, get standard deviation of distribution from diffusion constant
        stepStd = sqrt(2*diffCoeffitient*timeStep);

        % counter for timepoints
        countNumIterations = countNumIterations + numIterations;

        if constrained %if motion is confined within a certain region
        
         
          traj(countNumIterations-numIterations,2:dimension+1) = 0;
          
          for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points
                
                if i ==1
                    i = 2;
                end
    % for i=2:numIterations %go over all time points

        %get displacement
        displacementVec = stepStd*randn(1,dimension);

        %calculate square of new distance from center
        distance2 = sum((traj(i-1,2:dimension+1)+displacementVec).^2);

        %if new distance is larger than confinement radius, then bounce
        %the particle back
        if distance2 > cRadius^2
            
            prevPosVec = traj(i-1,2:dimension+1); %previoue position
            displacement2 = sum(displacementVec.^2); %square of mag. of original displacement

            %determine fraction of original displacement within confinement region
            dispFrac = roots([displacement2 2*displacementVec*prevPosVec' ...
                (sum(prevPosVec.^2)-cRadius^2)]); %solve equation
            dispFrac = dispFrac(find(dispFrac>0)); %assign correct root
            dispFrac = dispFrac(find(dispFrac<1)); %to fraction
            
            %get vector within confinement region and vector to be reflected
            vecInside = dispFrac*displacementVec;
            vecToReflect = (1-dispFrac)*displacementVec;
            
            %get radius vector and surface normal
            radiusVec = prevPosVec + vecInside;
            unitNormal = radiusVec/cRadius;

            %find component of vector (to be reflected) along surface normal
            normalComp = vecToReflect*unitNormal';

            %reflect vector
            vecToReflect = vecToReflect - 2*normalComp*unitNormal;
            
            %calculate new position of particle
            traj(i,2:dimension+1) = radiusVec + vecToReflect;
            
            %fix new position in case algorithm messed up due to some
            %numerical inaccuracies
            distance2 = sum(traj(i,2:dimension+1).^2);
            if distance2 > cRadius^2
                traj(i,2:dimension+1) = traj(i,2:dimension+1)*0.999*cRadius/sqrt(distance2);
            end

        else
            
            traj(i,2:dimension+1) = traj(i-1,2:dimension+1) + displacementVec;
        
        end %(if distance2 > cRadius^2)
        
    end %(for i=2:numIterations)
    
    if dimension == 1 || 2 || 3
     traj(countNumIterations-numIterations+1:countNumIterations,2) = traj(countNumIterations-numIterations+1:countNumIterations,2) + lastPosition(1);
     traj(countNumIterations-numIterations,2) = lastPosition(1);
    end
    
    if dimension == 2 || 3;
        
     traj(countNumIterations-numIterations+1:countNumIterations,3) = traj(countNumIterations-numIterations+1:countNumIterations,3) + lastPosition(2);
     traj(countNumIterations-numIterations,3) = lastPosition(2);
    end
    
    if dimension ==  3;
        
     traj(countNumIterations-numIterations+1:countNumIterations,4) = traj(countNumIterations-numIterations+1:countNumIterations,4) + lastPosition(3);
     traj(countNumIterations-numIterations,4) = lastPosition(3);
    end
    
     % saving last position of movement
     lastPosition = traj(i,2:dimension+1);
     
        elseif directed %if motion is directed

            for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points

                 if i ==1
                    i = 2;
                 end
                
                %get particle's position at this time point
                traj(i,2:dimension+1) = traj(i-1,2:dimension+1) + driftVelCart*timeStep ...
                    + stepStd*randn(1,dimension);

            end %if directed
            
        % saving last position of movement
        lastPosition = traj(i,2:dimension+1);

        else % if motion is simple diffusion

            for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points
                
                if i ==1
                    i = 2;
                end
                 
                traj(i,2:dimension+1) = traj(i-1,2:dimension+1)+ stepStd*randn(1,dimension);
            end %if simple diffusion
            
             
            
            
%              if dimension == 1 || 2 || 3
%      traj((countNumIterations-numIterations+1):countNumIterations,2) = traj((countNumIterations-numIterations+1):countNumIterations,2) + lastPosition(1);
%     end
%     
%     if dimension == 2 || 3;
%         
%      traj((countNumIterations-numIterations+1):countNumIterations,3) = traj((countNumIterations-numIterations+1):countNumIterations,3) + lastPosition(2);
%     end
%     
%     if dimension ==  3;
%         
%      traj((countNumIterations-numIterations+1):countNumIterations,4) = traj((countNumIterations-numIterations+1):countNumIterations,4) + lastPosition(3);
%     end
        % saving last position of movement
        lastPosition = traj(i,2:dimension+1);

        end %(if constrained)
end

%%%%% ~~ the end ~~ %%%%%

            
%             for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points
%                 
%                 if i ==1
%                     i = 2;
%                 end
%                 
%                 %get displacement
%                 displacementVec = stepStd*randn(1,dimension);
%                 
%                 %calculate square of new distance from center
%                 distance2 = sum((traj(i-1,2:dimension+1)+displacementVec).^2);
% 
%                 %if new distance is larger than confinement radius, then bounce
%                 %the particle back
%                 if distance2 > sum((cRadius + lastPosition).^2)
% 
%                     prevPosVec = traj(i-1,2:dimension+1)-lastPosition; %previoue position
%                     displacement2 = sum(displacementVec.^2); %square of mag. of original displacement
%                     
%                     
%                     
%                     %determine fraction of original displacement within confinement region
%                     dispFrac = roots([displacement2 2*displacementVec*prevPosVec' ...
%                         (sum(prevPosVec.^2)-cRadius^2)]); %solve equation
%                     dispFrac = dispFrac(find(dispFrac>0)); %assign correct root
%                     dispFrac = dispFrac(find(dispFrac<1)); %to fraction
% 
%                     %get vector within confinement region and vector to be reflected
%                     vecInside = dispFrac*displacementVec;
%                     vecToReflect = (1-dispFrac)*displacementVec;
% 
%                     %get radius vector and surface normal
%                     radiusVec = prevPosVec + vecInside;
%                     unitNormal = radiusVec/cRadius;
% 
%                     %find component of vector (to be reflected) along surface normal
%                     normalComp = vecToReflect*unitNormal';
% 
%                     %reflect vector
%                     vecToReflect = vecToReflect - 2*normalComp*unitNormal;
% 
%                     %calculate new position of particle
%                     traj(i,2:dimension+1) = radiusVec + vecToReflect +  traj(i,2:dimension+1);
% 
%                     %fix new position in case algorithm messed up due to some
%                     %numerical inaccuracies
%                     distance2 = sum(traj(i,2:dimension+1).^2);
%                     if distance2 > cRadius^2
%                         traj(i,2:dimension+1) = traj(i,2:dimension+1)*0.999*cRadius/sqrt(distance2);
%                     end
% 
%                 else
% 
%                     traj(i,2:dimension+1) = traj(i-1,2:dimension+1) + displacementVec;
% 
%                 end %(if distance2 > cRadius^2)
%                 
%             end %(for i=(countNumIterations-numIterations+1):countNumIterations)
%             
%             % saving last position of movement
%             lastPosition = traj(i,2:dimension+1);
%             
%         elseif directed %if motion is directed
% 
%             for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points
% 
%                 %get particle's position at this time point
%                 traj(i,2:dimension+1) = traj(i-1,2:dimension+1) + driftVelCart*timeStep ...
%                     + stepStd*randn(1,dimension);
% 
%             end %if directed
%             
%         % saving last position of movement
%         lastPosition = traj(i,2:dimension+1);
% 
%         else % if motion is simple diffusion
% 
%             for i=(countNumIterations-numIterations+1):countNumIterations %go over all time points
% 
%                 traj(i,2:dimension+1) = traj(i,2:dimension+1)+ stepStd*randn(1,dimension);
%             end %if simple diffusion
%             
%         % saving last position of movement
%         lastPosition = traj(i,2:dimension+1);
% 
%         end %(if constrained)
% end
% 
% %%%%% ~~ the end ~~ %%%%%
