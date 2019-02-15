function [diffModeDivider,fracTruePos] = estimDiffModeDividers(diffModeCoef,...
    numTraj,trajLength,doPlot,coordStd)
%ESTIMDIFFMODEDIVIDERS estimates the optimal dividers between different diffusion modes
%
%SYNOPSIS [diffModeDivider,fracTruePos] = estimDiffModeDividers(diffModeCoef,...
%    numTraj,trajLength,doPlot,coordStd)
%
%INPUT  diffModeCoef: Vector listing the diffusion coefficient of the
%                     different diffusion modes. The diffusion modes should
%                     be sorted in ascending order of their diffusion
%                     coefficients.
%       numTraj     : Number of trajectories to be simulated.
%                     Optional. Default: 1000.
%       trajLength  : Vector of simulated trajectory lengths, i.e. number
%                     of time points per trajectory.
%                     Optional. Default: [5:20].
%       doPlot      : 1 to plot the results, 0 otherwise.
%                     Optional. Default: 0.
%                     NO PLOTTING AT THE MOMENT.
%       coordStd    : Vector of standard deviations of coordinates (i.e.
%                     localization precision).
%                     Optional. Default: 0.
%
%OUTPUT diffModeDivider: Estimated dividers between the different modes
%                        (1st dimension) and their standard deviations (2nd
%                        dimension) for each coordinate standard deviation
%                        (3rd dimension) and each trajectory length (4th
%                        dimension).
%       fracTruePos    : Fraction of each mode correctly classified using
%                        the estimated dividers (1st dimension) and its
%                        standard deviation (2nd dimension) for each
%                        coordinate standard deviation (3rd dimension) and
%                        each trajectory length (4th dimension).
%
%REMARKS Code is written for 2D case only, but can be generalized to 3D.
%
%Khuloud Jaqaman, July 2012

%% Output
diffModeDivider = [];
fracTruePos = [];

%% Input

if nargin < 1
    disp('estimDiffModeDividers: Please enter vector of diffusion modes')
    return
end

if nargin < 2 || isempty(numTraj)
    numTraj = 1000;
end

if nargin < 3 || isempty(trajLength)
    trajLength = 5 : 20;
end

if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 5 || isempty(coordStd)
    coordStd = 0;
end

%% Divider estimation

%number of times to repeat estimation in order to get an std
numRep = 50;

%number of trajectory lengths
numLength = length(trajLength);

%number of coordinate standard deviation values
numCoordStd = length(coordStd);

%get number of diffusion modes
numMode = length(diffModeCoef);

%resere memory for output
diffModeDivider = NaN(numMode-1,2,numCoordStd,numLength);
fracTruePos = NaN(numMode-1,2,numCoordStd,numLength);

%go over all lengths
for iLength = 1 : numLength
    
    %get current trajectory length
    trajLengthC = trajLength(iLength);
    
    %reserve memory for intermediate results
    diffModeDividerRep = NaN(numMode-1,numRep,numCoordStd);
    fracTruePosRep = NaN(numMode-1,numRep,numCoordStd);
    
    %go over repetitions
    progressText(0,'Estimating dividers');
    for iRep = 1 : numRep
        
        %simulate trajectories in each mode
        traj0 = NaN(trajLengthC,2,numTraj,numMode);
        for iMode = 1 : numMode
            for iTraj = 1 : numTraj
                trajTmp = brownianMotion(2,diffModeCoef(iMode),trajLengthC,0.01);
                trajTmp = trajTmp(1:100:end,:);
                trajTmp = trajTmp(2:end,:);
                traj0(:,:,iTraj,iMode) = trajTmp;
            end
        end
        
        for iCoordStd = 1 : numCoordStd
            
            %perturb positions to reflect localization error
            traj = traj0 + randn(size(traj0))*coordStd(iCoordStd);
            
            %calculate diffusion coefficient
            msd = squeeze( mean( sum( diff(traj).^2,2 ) ) );
            diffCoef = msd / 4 - coordStd(iCoordStd)^2;
            
            %for each pair of consecutive modes, slide divider between their means and
            %calculate the fraction of each mode that will be properly classified
            %given divider value
            divValue = 0:0.001:diffModeCoef(end)*1.1;
            numValue = length(divValue);
            fracModeBelow = NaN(numValue,numMode);
            fracModeAbove = NaN(numValue,numMode);
            for iMode = 1 : numMode
                for iValue = 1 : numValue
                    fracModeBelow(iValue,iMode) = length(find(diffCoef(:,iMode)<=divValue(iValue)))/numTraj;
                    fracModeAbove(iValue,iMode) = length(find(diffCoef(:,iMode)>divValue(iValue)))/numTraj;
                end
            end
            
            %find the best divider between each pair of modes, i.e. the divider value
            %which simultaneously maximizes the correct classification of each mode
            for iMode = 1 : numMode-1
                fracDiff = abs( fracModeBelow(:,iMode) - fracModeAbove(:,iMode+1) );
                dividerIndx = find(fracDiff == min(fracDiff));
                dividerIndx = round(median(dividerIndx));
                diffModeDividerRep(iMode,iRep,iCoordStd) = divValue(dividerIndx);
                fracTruePosRep(iMode,iRep,iCoordStd) = fracModeBelow(dividerIndx,iMode);
            end
            
        end
        
        progressText(iRep/numRep,'Estimating dividers');
        
    end %(for iRep = 1 : numRep)
    
    %calculate the average and std for output
    diffModeDivider(:,1,:,iLength) = mean(diffModeDividerRep,2);
    diffModeDivider(:,2,:,iLength) = std(diffModeDividerRep,[],2);
    fracTruePos(:,1,:,iLength) = mean(fracTruePosRep,2);
    fracTruePos(:,2,:,iLength) = std(fracTruePosRep,[],2);
    
end %(for iLength = 1 : numLength)

%% Plotting

% %plot the fractions and display dividers
% if doPlot
%     
%     %for display, remove unnecessary regions
%     fracModeBelowDisp = fracModeBelow;
%     fracModeAboveDisp = fracModeAbove;
%     for iMode = 1 : numMode
%         fracModeBelowDisp(divValue<diffModeCoef(iMode),iMode) = NaN;
%         fracModeAboveDisp(divValue>diffModeCoef(iMode),iMode) = NaN;
%         indx = find(fracModeBelowDisp(:,iMode)==1,6,'first');
%         if ~isempty(indx)
%             fracModeBelowDisp(indx(end):end,iMode) = NaN;
%         end
%         indx = find(fracModeAboveDisp(:,iMode)==1,6,'last');
%         if ~isempty(indx)
%             fracModeAboveDisp(1:indx(1),iMode) = NaN;
%         end
%     end
%     
%     %display the dividers on the last simulation results
%     figure, hold on
%     plot(divValue,fracModeBelowDisp(:,1:end-1))
%     plot(divValue,fracModeAboveDisp(:,2:end))
%     plot(diffModeDivider(:,1),fracTruePos(:,1),'ko')
%     
% end

%% ~~~ the end ~~~


% %display a histogram of the diffuion coefficients and indicate mode
% %dividers
% figure
% hist(log10(diffCoef),100)
% hold on
% for iMode = 1 : numMode-1
%     plot(log10(diffModeDivider(iMode))*[1 1],[0 250])
% end
%
% %calculate fraction of each mode included within dividers
% modeFracAboveBelowDivider = NaN(numMode,2);
% modeFracAboveBelowDivider(1,2) = length(find(diffCoef(:,1)<diffModeDivider(1)))/size(diffCoef,1);
% for i = 2 : numMode-1
%     modeFracAboveBelowDivider(i,1) = length(find(diffCoef(:,i)>diffModeDivider(i-1)))/size(diffCoef,1);
%     modeFracAboveBelowDivider(i,2) = length(find(diffCoef(:,i)<diffModeDivider(i)))/size(diffCoef,1);
% end
% i = numMode;
% modeFracAboveBelowDivider(i,1) = length(find(diffCoef(:,i)>diffModeDivider(i-1)))/size(diffCoef,1);
