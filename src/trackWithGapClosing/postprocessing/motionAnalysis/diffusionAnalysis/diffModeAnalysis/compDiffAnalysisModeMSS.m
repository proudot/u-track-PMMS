function [compMode2MSS,compModeWithMSS2ModeNoMSS] = compDiffAnalysisModeMSS(tracksFinal,...
    diffAnalysisMSS,diffAnalysisMode,lengthMinMax)
%COMPDIFFANALYSISMODEMSS compares mode decomposition to MSS diffusion analysis
%
%SYNOPSIS compDiffAnalysisModeMSS(tracksFinal,diffAnalysisMSS,diffAnalysisMode)
%
%INPUT  tracksFinal     : Output of trackCloseGapsKalman.
%       diffAnalysisMSS : Output of trackDiffusionAnalysis1.
%       diffAnalysisMode: Output of trackDiffModeAnalysis.
%       lengthMinMax    : Minimum and maximum length of trajectories to
%                         include in analysis.
%                         Optional. Default: [5 99]
%
%OUTPUT 
%
%Khuloud Jaqaman, June 2013

%% Input

if nargin < 3
    disp('--compDiffAnalysisModeMSS: Missing input arguments!');
    return
end

if nargin < 4 || isempty(lengthMinMax)
    lengthMinMax = [5 99];
end


%% Pre-processing

%keep only tracks of appropriate length
criteria.lifeTime.min = lengthMinMax(1);
criteria.lifeTime.max = lengthMinMax(2);
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx);
diffAnalysisMSS = diffAnalysisMSS(indx);
diffAnalysisMode = diffAnalysisMode(indx);

%extract diffusion analysis results
diffClassMSS = vertcat(diffAnalysisMSS.classification);
diffClassMSS = diffClassMSS(:,2);
diffClassMode = vertcat(diffAnalysisMode.diffMode);
diffCoefMSS = catStruct(1,'diffAnalysisMSS.fullDim.genDiffCoef(:,3)');
diffCoefMode = vertcat(diffAnalysisMode.diffCoef);

%get indices of tracks with MSS analysis and Mode analysis, as well as
%indices of tracks with only Mode analysis (i.e. no MSS analysis)
indxMSS  = find( ~isnan(diffClassMode) & ~isnan(diffClassMSS) );
indxMode = find( ~isnan(diffClassMode) & isnan(diffClassMSS) );
numTrajMSS = length(indxMSS);
numTrajMode = length(indxMode);

%keep only the information of those tracks
diffClassMSS_MSS = diffClassMSS(indxMSS);
diffClassMode_MSS = diffClassMode(indxMSS);
diffCoefMSS_MSS = diffCoefMSS(indxMSS);
diffCoefMode_MSS = diffCoefMode(indxMSS);
diffClassMode_Mode = diffClassMode(indxMode);
diffCoefMode_Mode = diffCoefMode(indxMode);

%get numbers
numMode = max(diffClassMode);
numMSS = 3;
numTracks = length(tracksFinal);

%% Classification correspondence and probabilities

%calcualte absolute probabilities, i.e. sum(all probabilities) = 1
probModeMSSAbs = zeros(numMode,numMSS);
for iMode = 1 : numMode
    for iMSS = 1 : numMSS
        probModeMSSAbs(iMode,iMSS) = length( find( diffClassMode_MSS==iMode & ...
            diffClassMSS_MSS==iMSS ) ) / numTrajMSS;
    end
end

%calculate relative probebailities within each mode, i.e. sum(probabilities
%within each mode) = 1
probModeMSSRelPerMode = probModeMSSAbs ./ repmat(sum(probModeMSSAbs,2),1,numMSS);
probModeMSSRelPerMSS = probModeMSSAbs ./ repmat(sum(probModeMSSAbs,1),numMode,1);

%compare the distribution of modes for tracks long enough for MSS analysis
%to the distribution of modes for tracks that are too short for MSS analysis
probModeWithMSS = hist(diffClassMode_MSS,1:numMode);
probModeWithMSS = probModeWithMSS / sum(probModeWithMSS);
probModeNoMSS = hist(diffClassMode_Mode,1:numMode);
probModeNoMSS = probModeNoMSS / sum(probModeNoMSS);

%% Diffusion coefficient correspondence

%get the positional error in each track
coordVar = [];
for iTrack = 1 : 1000 : numTracks
    trackMat = convStruct2MatIgnoreMS(tracksFinal(iTrack:min(iTrack+1000,numTracks)),1);
    xCoordStd = trackMat(:,5:8:end);
    yCoordStd = trackMat(:,6:8:end);
    coordVar = [coordVar; nanmean([xCoordStd yCoordStd].^2,2)];
end
meanPosVar = nanmean(coordVar);

%fit a straight line through the Mode diff coef vs. MSS diff coef, after
%correcting the MSS diff coef for positional error
[lineCoefVal,S] = polyfit(diffCoefMode_MSS,diffCoefMSS_MSS-meanPosVar,1);
Rinv = inv(S.R);
varCovMat = (Rinv*Rinv')*(S.normr^2)/S.df;
lineCoefStd = sqrt(diag(varCovMat));

%compare the distribution of diffusion coefficients for tracks long enough
%for MSS analysis to teh distribution for tracks that are too short for MSS
%analysis
pvKS = zeros(500,1);
for i=1:500
    indx1 = randsample(numTrajMSS,400);
    indx2 = randsample(numTrajMode,400);
    [~,pvKS(i)] = kstest2(diffCoefMode_MSS(indx1),diffCoefMode_Mode(indx2));
end
pvKS = mean(pvKS);

%calculate the mean and std of the diffusion coefficient distributions
diffCoefModeWithMSSMeanStd = [mean(diffCoefMode_MSS) std(diffCoefMode_MSS)];
diffCoefModeNoMSSMeanStd = [mean(diffCoefMode_Mode) std(diffCoefMode_Mode)];

%% Output
compMode2MSS = struct('probModeMSSAbs',probModeMSSAbs,...
    'probModeMSSRelPerMode',probModeMSSRelPerMode,...
    'probModeMSSRelPerMSS',probModeMSSRelPerMSS,...
    'diffCoefRegLine',[lineCoefVal' lineCoefStd]);
compModeWithMSS2ModeNoMSS = struct('probMode',[probModeWithMSS;probModeNoMSS],...
    'diffCoefMeanStd',[diffCoefModeWithMSSMeanStd;diffCoefModeNoMSSMeanStd],...
    'pvDiffCoefKS',pvKS);

%% ~~~ the end ~~~

