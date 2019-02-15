function [trackBehavior] = transBehaviorAnalysis...
         (trajectoryMatrix, tWinSize, nDiff, nDev, nMSD)
     
% TRANSBEHAVIORANALYSIS analysis the transient behavior of trajectories by
%                       calculating the diffusion coeffitient, the
%                       deviation of the MSD and the assymetry with a 
%                       rolling time window of varying size for the
%                       different parameters.
%
% SYNOPSIS [trackBehavior] = transBehaviorAnalysis...
%          (trajectoryMatrix, tWinSize, nDiff, nDev, nMSD)
% 
% INPUT  trajectoryMatrix = n-by-8 matrix with n trajectories and columns
%                            (x, y, intensity, amplitude,
%                             dx, dx, dintensity, damplitude)
%
%        tWinSize         = 3-by-2 matrix with 
%                             rows: 1. for diffusion
%                                   2. for constrained
%                                   3. for directed
%                             columns : 1. minimum Size of the rolling time
%                                          window
%                                       2. maximum size of the rolling time
%                                          window       
%                           Guidelines for choosing the time window size:
%                         - the minimum size of the time window should be
%                           smaller then the minimum duration of a certain
%                           behavior. This can only be rigorously verified
%                           a posteriori to the analysis.
%                         - the mean duration of a trajectory could be used
%                           as the maximum window size. It can be reduced
%                           by a priori estimation with a test pool of 
%                           trajectoriesto save computational time
%
%        nDiff (optional) = number of timpoints of MSD which are used for
%                           the linear regression fit to determine the 
%                           diffussion coeffitions. (default = 5)
%
%        nDev (optional)  = number of timepoints of MSD compared to the
%                           linear regression fit to determine the
%                           deviation between them. (default = 50)
%
%        nMSD (optional)  = number of timepoints used for the calculation
%                           of the MSD (min = 4, max = tWinSizeMin-2,
%                           default = 20) 
%
% OUTPUT  trackBehavior   = n-by-18-by-l matrix with
%             n rows:  timepoints of trajectories
%             columns: 1. minimum diffusion coeffitient of all different 
%                         time windows width for each time point
%                      2. p-Value, which reflects the probability that
%                         sigma Hat Square is fischer distributed
%                         with a mean of 1.
%                      3. minimum deviation of the MSD from linearity (fit 
%                         over first nDiff points of MSD) of all different 
%                         time windows width for each time point
%                      4. maximum asymmetry of all different time windows 
%                         width for each time point 
%                    5-7. empty (NaN)
%                      8. window size used for computing the min diffusion
%                         coeffitient in column 1
%                      9. window size used for computing the min deviation
%                         of the MSD from linearity in column 3
%                     10. window size used for computing the max asymmetry
%                         in column 3
%                  11-15. empty (NaN)
%                     16. actual number of data points used for the
%                         calculation of the min diffusion coeffitient
%                         (column 1)
%                     17. actual number of data points used for the
%                         calculation of the min deviation of the MSD from
%                         linearity (column 3)
%                     18. actual number of data points used for the
%                         calculation of the max asymmetry (column 4)
%             l: trajectories
%
% CREATED  gp 2/13/07

% ----------------------------------------------
% checking the input variables
% ----------------------------------------------

% checks if the number of time points to calculate the MSD (nMSD) is given
% by the user, if not it is set to default(20)
if (nargin <= 4);
    nMSD = 20;
end

% checks if the number of time points to calculate the deviation (nDev) is
% given by the user, if not it is set to default(5)
if (nargin <= 3);
    nDev = 50;
end

% checks if the number of time points to calculate the diffusion (nDiff) is
% given by the user, if not it is set to default(5)
if (nargin <= 2);
    nDiff = 5;
end

% checks if the time window sizes are in the right format
if sum(any(mod(tWinSize, 2) ~= 1));
    error('the minimun and the maximum time window sizes must be odd integers');
end

% checks that the minimun size of the rolling time window is bigger then 1
if any(tWinSize(:,1) <= 1);
    error('the minimum size of the rolling time window (tWinSizeMin) has to be >= 3');
end 

% determines the row and column size of the trajectory matrix
[rowSize,columnSize] = size(trajectoryMatrix);

% checks if the input variable tWinSizeMax is smaller than the number of
% timepoints in a given trajectory
if any(tWinSize(:,2)) >= (columnSize/8);
     error(['the maximum size of the rolling time window (tWinSizeMax) has' ...
            'to be smaller than the number of timepoints of a given trajectory']);
end

% checks if the number of timepoints in the smallest time window for the 
% calculation of the diffusion coeffitient is bigger then the number of time
% points needed by the function which calculates the MSD
if nMSD <= 4 || nMSD > (tWinSize(1,1)-2);
     error(['the number of timepoints for the MSD calculation  (nMSD) has',...
            'to be >= 4 and <= (tWinSizeDiffMin - 2). Please decrease nMSD or',...
            'if you want to keep the default settings, please increase',...
            'tWinSizeDiffMin over 21']);
end 

% checks if the input trajctory matrix is in the right format
if mod(columnSize,8) ~= 0;
     error(['the column size must be a multiplicate of 8 in the form of',...
            '[X-value(1), Y-value(2),dx(5), dy(6)]']);
end

%checks if the number of timepoints for the MSD calculation (nMSD) is
%increased by nDiff to ensure that enough MSDs are calculated for the
%diffusion calculation
if tWinSize(1,1) < nMSD + nDiff
    error('the minimum widow size has to be smaller then nDiff and nMSD together');
end

% ----------------------------------------------
% initializing
% ----------------------------------------------

% initialize sets output matrix to NaN;
trackBehavior= ones(columnSize/8,18,rowSize) *NaN;

% ----------------------------------------------
% assignment of variables
% ----------------------------------------------

% determines the minimum and the maximum size of all window sizes for the
% different parameters
tWinSizeMin = min(tWinSize(:,1));
tWinSizeMax = max(tWinSize(:,2));

% creates the input structure for the different parameters
analyzingParameters = struct('tWinSizeDiffMin',tWinSize(1,1),...
    'tWinSizeDiffMax',tWinSize(1,2), 'tWinSizeConMin',tWinSize(2,1),...
    'tWinSizeConMax',tWinSize(2,2),'tWinSizeDirMin',tWinSize(3,1),...
    'tWinSizeDirMax',tWinSize(3,2),'tWinSizeMin',tWinSizeMin,...
    'tWinSizeMax',tWinSizeMax,'nDiff', nDiff,'nDev',nDev, 'nMSD',nMSD);

% ----------------------------------------------
% computation of the different trajectories
% ----------------------------------------------

%determines the beginning and the end of data points in the trajectories
trackSEL = getTrackSEL(trajectoryMatrix);

% loops through the different trajectories (rows) of the given 
% trajectory matrix (trajectoryMatrix)  
for n=1:rowSize     
    
    % creates a matrix of one trajectory (row) and saves it into the 
    % variable trajectoryTimepointMatrix                                                                               
    trajectoryTimepointMatrix = trajectoryMatrix(n,trackSEL(n,1)*8-7:trackSEL(n,2)*8);  
    
    % reshapes the vector matrix (trajectoryTimepointMatrix) into a 
    % timepoint = row and value = column matrix                                                                                
    trajectoryTimepointMatrix = reshape(trajectoryTimepointMatrix,8,trackSEL(n,3))';  
    
    % checks if the trajectoryTimepointMatrix contains as least as many data
    % points as the smallest minimum window size for one of the three parameters
    if size(trajectoryTimepointMatrix,1) >= tWinSizeMin
        % calls the function transBehaviorTrack which computes the minimum 
        % diffusion coeffitient (first column), the corresponding p-Value 
        % (second column)of the sigma hat 
        % square (if it is fischer distributed with a mean of 1) and the 
        % minimum deviation between the linear fit of the MSD and the MSD of 
        % all given time windows for ech time point of a given trajectory.                                                                                                                                                               
        [trackBehavior(trackSEL(n,1):trackSEL(n,2),:,n)] =...
            transBehaviorTrack(trajectoryTimepointMatrix, analyzingParameters); 
    end
    
    % displays the number of calculated trajectories
    display(n);
end


