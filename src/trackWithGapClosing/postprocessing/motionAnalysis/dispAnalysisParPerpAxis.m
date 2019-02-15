function [meanDispPar,meanDispPerp,meanSqDispPar,meanSqDispPerp] = ...
    dispAnalysisParPerpAxis(tracks,minTrackLen,diffAnalysisRes,maxLag)
%dispAnalysisParPerpAxis analyzes the mean displacement and square displacement parallel and perpendicular to the direction of motion
%
%SYNOPSIS [meanDispPar,meanDispPerp,meanSqDispPar,meanSqDispPerp] = ...
%    dispAnalysisParPerpAxis(tracks,minTrackLen,diffAnalysisRes,maxLag)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in analysis.
%                    Optional. Default: 5.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffusionAnalysis1). Optional. If not input, will
%                    be calculated here.
%
%OUTPUT meanDispPar   : maxLag-by-3-by-4 array. Rows indicate mean
%                       parallel displacement at lags 1, 2, etc. Columns
%                       indicate mean, std and sem. Third dimension
%                       indicates track type (see below).
%       meanDispPerp  : maxLag-by-3-by-4 array. Rows indicate mean
%                       perpendicular displacement at lags 1, 2, etc. 
%                       Columns indicate mean, std and sem. Third dimension
%                       indicates track type (see below).
%       meanSqDispPar : maxLag-by-3-by-4 array. Rows indicate mean
%                       square parallel displacement at lags 1, 2, etc. 
%                       Columns indicate mean, std and sem. Third dimension
%                       indicates track type (see below).
%       meanSqDispPerp: maxLag-by-3-by-4 array. Rows indicate mean
%                       square perpendicular displacement at lags 1, 2, etc. 
%                       Columns indicate mean, std and sem. Third dimension
%                       indicates track type (see below).
%
%Khuloud Jaqaman, October 2010

%% output
meanDispPar = NaN(maxLag,3,4);
meanDispPerp = NaN(maxLag,3,4);
meanSqDispPar = NaN(maxLag,3,4);
meanSqDispPerp = NaN(maxLag,3,4);

%% input

if nargin < 2 || isempty(tracks)
    disp('dispAnalysisParPerpAxis: Missing input arguments!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(diffAnalysisRes)
    diffAnalysisRes = trackDiffusionAnalysis1(tracks,1,2,1,[0.05 0.2]);
end

if nargin < 4 || isempty(maxLag)
    maxLag = 10;
end

%% preamble

%ignore merges and splits and divide compound tracks back into the
%individual track segments
inputStruct = tracks;
clear tracksFinal
tracks = convStruct2MatIgnoreMS(inputStruct);

%extract track classification from diffAnalysisRes
trackType = vertcat(diffAnalysisRes.classification);

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
tracks = tracks(indx,:);
trackType = trackType(indx,:);

%find linear tracks with normal diffusive behavior along preferred
%direction of motion
tracks1 = tracks(trackType(:,1)==1 & trackType(:,3)==2,:);

%find linear tracks with super-diffusive behavior along preferred direction
%of motion
tracks2 = tracks(trackType(:,1)==1 & trackType(:,3)==3,:);

%find non-linear tracks with normal diffusive behavior
tracks3 = tracks(trackType(:,1)==0 & trackType(:,2)==2,:);

%find non-linear tracks with sub-diffusive behavior
tracks4 = tracks(trackType(:,1)==0 & trackType(:,2)==1,:);

clear criteria tracks indx

%% displacement statistics

%go over all types of linear tracks
for iType = 1 : 4
    
    %get information on current track type and reserve memory
    eval(['tracks = tracks' num2str(iType) ';']);
    numTracks = size(tracks,1);
    
    %go over all lags
    for iLag = 1 : maxLag
        
        %initialize variables
        dispParIn = [];
        dispPerpIn = [];
        
        %go over all tracks in this category
        for iTrack = 1 : numTracks
            
            %get the positions in this track and their standard deviations
            %keep NaNs to mark gaps
            trackCoordX = tracks(iTrack,1:8:end)';
            deltaCoordX = tracks(iTrack,5:8:end)';
            trackCoordY = tracks(iTrack,2:8:end)';
            deltaCoordY = tracks(iTrack,6:8:end)';
            trackCoord = [trackCoordX trackCoordY];
            deltaCoord = [deltaCoordX deltaCoordY];
            
            %project positions onto track's direction of motion
            [~,~,velDir] = projectCoordOntoDir(...
                trackCoord,deltaCoord,[],[]);
            
            %calculate vector of displacements
            %keep NaNs to mark gaps
            dispVec = trackCoord(1+iLag:end,:) - trackCoord(1:end-iLag,:);
            
            %calculate the dot product of displacements with the direction
            %vector and its normal
            dispAlongDir = dispVec * velDir;
            dispPerpDir = dispVec * [velDir(2) -velDir(1)]';
            
            %collect all displacement magnitudes parallel to direction of
            %motion
            dispParT = abs(dispAlongDir(~isnan(dispAlongDir)));
            
            %collect all displacement magnitudes perpendicular to direction of
            %motion
            dispPerpT = abs(dispPerpDir(~isnan(dispPerpDir)));
            
            %add this track's results to the rest
            dispParIn = [dispParIn; dispParT];
            dispPerpIn = [dispPerpIn; dispPerpT];
            
        end %(for iTrack = 1 : numTracks)
        
        meanDispPar(iLag,:,iType)  = [mean(dispParIn) std(dispParIn) length(dispParIn)];
        meanDispPerp(iLag,:,iType) = [mean(dispPerpIn) std(dispPerpIn) length(dispPerpIn)];
        meanSqDispPar(iLag,:,iType)  = [mean(dispParIn.^2) std(dispParIn.^2) length(dispParIn)];
        meanSqDispPerp(iLag,:,iType) = [mean(dispPerpIn.^2) std(dispPerpIn.^2) length(dispPerpIn)];
        
    end
    
end %(for iType = 1 : 4)

%% ~~~ the end ~~~



