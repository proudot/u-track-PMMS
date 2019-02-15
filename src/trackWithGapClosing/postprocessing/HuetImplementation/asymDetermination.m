function asymmetry = asymDetermination(positions);

% ASYMDETERMINATION  function which calculates the asymmetry of the gyration
%                    tensor of the x and y values of a given set of
%                    datapoints (trajectory)
%
% SYNOPSIS  asymmetry = asymDetermination(positions)
% 
% INPUT     positions = n-by-m matrix
%                   
%                   rows: timepoints along a trajectory
%                   columns: 1. x-values
%                            2. y-values  
%
% OUTPUT    asymmetry = asymmetry of all positions in the dataset
%                       (trajectory)
%
% CREATED gp 4/02/07

%-----------------------------------------
% Initialize
%-----------------------------------------

asymmetry = [];

% find ~NaNs and store the location in goodData. This is done to supress
% errors which will occur if NaNs are computed for asymmetry
goodData = find(isfinite(positions(:,1)) & isfinite(positions(:,2)));

% executes the asymmetry calculation if at least 2/3 of the datapoints
% are ~NaN
if size(goodData,1) >= 2/3 * size(positions,1);

    %----------------------------------------------------------
    % Calculation of the 2D-radius (Rg) of the gyration tensor 
    %----------------------------------------------------------

    % creating position 1,1 of the Rg matrix
    R11 = (mean(positions(goodData,1).^2)) - ((mean(positions(goodData,1))).^2);

    % creating position 2,2 of the Rg matrix
    R22 = (mean(positions(goodData,2).^2)) - ((mean(positions(goodData,2))).^2);

    % creating position 1,2 of the Rg matrix
    R12 = (mean(positions(goodData,1).*positions(goodData,2))) - (mean(positions(goodData,1)) * (mean(positions(goodData,2))));

    % setting position 2,1 as position 1,2 of the Rg matrix
    R21 = R12;

    % creation of the Rg matrix
    Rg = [R11 R12; R21 R22];

    % computes the eigen value and the eigenvector of gyration tensor Rg
    [eigenVector, eigenValue] = eig(Rg);

    % computes the gyration Radius for the first dimension Rxx 
    Rxx = sqrt(eigenValue(1,1));

    % computes the gyration Radius for the second dimension Ryy 
    Ryy = sqrt(eigenValue(2,2));

    % calculation of the asymmetry 
    asymmetry = -log(1-(((Rxx^2 - Ryy^2)^2) / (Rxx^2 + Ryy^2)^2));
    
else % if there are more then 1/3 of NaNs, the asymmetry will not be 
     %calculated and gives the output NaN
    
    asymmetry = NaN;
    
end % of if size(goodData,1) >= 2/3 * size(positions,1);