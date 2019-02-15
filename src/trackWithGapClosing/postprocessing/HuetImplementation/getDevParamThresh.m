nDiff = 5;
nMSD = 20;
nDev = 50;
numSamples = 1000;
numTimePoints = 501;

%traj1 = NaN*ones(numTimePoints,2,numSamples);
traj0p1 = NaN*ones(numTimePoints,2,numSamples);
%Dev1 = NaN*ones(numTimePoints,numSamples);
Dev0p1 = NaN*ones(numTimePoints,numSamples);

%generate numSample trajectories
for iSample = 1 : numSamples

%     %generate 100 time point random walk with diffusion coefficient = 1 micron^2/s
%     [traj1(:,:,iSample),errFlag] = brownianMotion(2,1,50,0.1,1,1);

    %generate 100 time point random walk with diffusion coefficient = 0.1 micron^2/s
    [traj0p1(:,:,iSample),errFlag] = brownianMotion(2,0.5,50,0.1);

end


% %calculate the deviation parameter for the generated trajectories
% %go from length 3 time points to total length in the 1 micron^2/s
% %trajectory
% for iSample = 1 : numSamples
%     for iTP = 501 : 50 : numTimePoints
% %        Asym1(iTP,iSample) = asymDetermination(traj1(1:iTP,:,iSample));
%        
%         %----------------------------------------------------
%         % Convertion of the timeWinMatrix
%         %----------------------------------------------------
%         % converts the timeWinMatrix into a structure array which is the
%         % input format for the meanSquaredDisplacement function.
% 
%         % coordinates include the x and y values of the timeWinMatrix 
%         positions.coordinates = [traj1(1:iTP,1,iSample),traj1(1:iTP,2,iSample)];
% 
%         % creates the 2-by-2-by-rownumber stucture matrix 
%         % positions.covariances and sets all values to zero. 
%         positions.covariances = zeros(2,2,iTP);
% 
%         % fills the structure matrix positions.covariances with the sigmaX2
%         % values from the timeWinMatrix
%         positions.covariances(1,1,:) = 0.000001;
% 
%         % fills the structure matrix positions.covariances with the sigmaY2
%         % values from the timeWinMatrix
%         positions.covariances(2,2,:) = 0.000001;
% 
%         % creates the input structure for the different parameters
%         analyzingParameters = struct('nDiff', nDiff,'nMSD',nMSD, 'nDev', nDev);
% 
%         % calles the function diffussion which computes the 
%         % diffussion coeffitient per time point over the given window size
%         [diffCoeffitientForDev, pValueForDev,mSqDisp, linSlope,...
%             yIntercept] = diffussion(positions,analyzingParameters);   
% 
%         % calles the function standardDev which calculates the mean deviation
%         % between the regression fit of the MSD (*nDiff) and the MSD of the 
%         % different timepoints. 
%         Dev1(iTP,iSample) = linearDev(analyzingParameters,mSqDisp,linSlope, yIntercept);
%      
%     end
%     display(iSample);
% end

% calculate the asymmetry parameter for the generated trajectories
% go from length 3 time points to total length in the 0.1 micron^2/s
% trajectory
for iSample = 1 : numSamples
    for iTP = 75 :150:numTimePoints
%        asym0p1(iTP,iSample) = asymDetermination(traj0p1(1:iTP,:,iSample));

        %----------------------------------------------------
        % Convertion of the timeWinMatrix
        %----------------------------------------------------
        % converts the timeWinMatrix into a structure array which is the
        % input format for the meanSquaredDisplacement function.

        % coordinates include the x and y values of the timeWinMatrix 
        positions.coordinates = [traj0p1(1:iTP,1,iSample),traj0p1(1:iTP,2,iSample)];

        % creates the 2-by-2-by-rownumber stucture matrix 
        % positions.covariances and sets all values to zero. 
        positions.covariances = zeros(2,2,iTP);

        % fills the structure matrix positions.covariances with the sigmaX2
        % values from the timeWinMatrix
        positions.covariances(1,1,:) = 0.000001;

        % fills the structure matrix positions.covariances with the sigmaY2
        % values from the timeWinMatrix
        positions.covariances(2,2,:) = 0.000001;

        % creates the input structure for the different parameters
        analyzingParameters = struct('nDiff', nDiff,'nMSD',nMSD, 'nDev', nDev);

        % calles the function diffussion which computes the 
        % diffussion coeffitient per time point over the given window size
        [diffCoeffitientForDev, pValueForDev,mSqDisp, linSlope,...
            yIntercept] = diffussion(positions,analyzingParameters);   

        % calles the function standardDev which calculates the mean deviation
        % between the regression fit of the MSD (*nDiff) and the MSD of the 
        % different timepoints. 
        Dev0p1(iTP,iSample) = linearDev(analyzingParameters,mSqDisp,linSlope, yIntercept);    
        
        % calculation of devRel which determines the relative deviation
        % corresponding to the min window size which was used to calculate the
        % minimum deviation of each time point 
        Dev99(iTP,iSample) = Dev0p1(iTP,iSample) - ((-0.0097*(log(iTP).^2))...
            + 0.1413*(log(iTP)) -1.1085);

    end
    display(iSample);

end

% %get the asymmetry parameter threshold, such that 99% of the sample has a
% %lower assymmetry parameter
% thresh1 = prctile(asym1,99,2);
% thresh0p1 = prctile(asym0p1,99,2);

%get the deviation threshold, such that 99% of the sample has a
%lower assymmetry parameter
%thresh1 = prctile(Dev1,1,2);
thresh99 = prctile(Dev99,7,2);