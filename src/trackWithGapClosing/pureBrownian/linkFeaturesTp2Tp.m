function [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
    linkFeaturesTp2Tp(movieInfo,costMatFun,costMatParams)
%LINKFEATURESTP2TP links features between consecutive time points in a movie using LAP
%
%SYNOPSIS [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
%    linkFeaturesTp2Tp(movieInfo,costMatFun,costMatParams)
%
%INPUT  movieInfo    : Array of size equal to the number of time points
%                      in a movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%             .yCoord      : Image coordinate system y-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D. Default: zeros.
%             .zCoord      : Image coordinate system z-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D or 2D. Default: zeros.
%             .amp         : Amplitudes of PSFs fitting detected features. 
%                            1st column for values and 2nd column 
%                            for standard deviations.
%       costMatFun   : Name of function used to calculate the cost matrix
%                      for linking. 'costMatSimple' for simple cost matrix,
%                      'costMatLogL' for cost matrix that uses information
%                      on the distribution of diplacements and amplitudes.
%       costMatParams: Structure with fields providing the parameters needed 
%                      to calculate the cost matrix. See cost matrix of 
%                      interest to determine what fields to include.
%
%OUTPUT trackedFeatureNum: Connectivity matrix of features between time points.
%                          Rows indicate continuous tracks, while columns 
%                          indicate time points. A track that ends before the
%                          last time point is followed by zeros, and a track
%                          that starts at a time after the first time point
%                          is preceded by zeros. 
%       trackedFeatureInfo:The positions and amplitudes of the tracked
%                          features. Number of rows = number of tracks, 
%                          while number of columns = 8*number of time 
%                          points. Each row consists of 
%                          [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                          in image coordinate system (coordinates in
%                          pixels). NaN is used to indicate time points 
%                          where the track does not exist.
%       errFlag          : 0 if function executes normally, 1 otherwise.
%
%REMARKS No gap closing.
%        The algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must have some
%        features in it.
%
%Khuloud Jaqaman, August 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackedFeatureNum = [];
trackedFeatureInfo = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('linkFeaturesTp2Tp')
    disp('--linkFeaturesTp2Tp: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get number of time points in movie
numTimePoints = length(movieInfo);

%check whether problem is 1D, 2D or 3D and augment coordinates if necessary
if ~isfield(movieInfo,'yCoord') %if y-coordinates are not supplied

    %assign zeros to y and z coordinates
    for i=1:numTimePoints
        movieInfo(i).yCoord = zeros(size(movieInfo(i).xCoord));
        movieInfo(i).zCoord = zeros(size(movieInfo(i).xCoord));
    end

else %if y-coordinates are supplied

    if ~isfield(movieInfo,'zCoord') %if z-coordinates are not supplied

        %assign zeros to z coordinates
        for i=1:numTimePoints
            movieInfo(i).zCoord = zeros(size(movieInfo(i).xCoord));
        end

    end %(if ~isfield(movieInfo,'zCoord'))

end %(if ~isfield(movieInfo,'yCoord') ... else ...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Linking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of features at each time point
for t=1:numTimePoints
    movieInfo(t).num = size(movieInfo(t).xCoord,1);
end

%fill the feature numbers in first time point in the connectivity matrix
trackedFeatureNum = [1:movieInfo(1).num]';

%go over all time points
for t = 1:numTimePoints-1

    %get the number of features in the two time points
    n = movieInfo(t).num;
    m = movieInfo(t+1).num;

    if n ~= 0 %if there are features in first frame

        if m ~= 0 %if there are features in second frame

            %calculate cost matrix
            eval(['[costMat,noLinkCost,nonlinkMarker] = ' costMatFun ...
                '(movieInfo(t:t+1),costMatParams);']);

            if any(costMat(:)~=nonlinkMarker) %if there are potential links
            
                %track features based on this cost matrix, allowing for birth and death
                [link12,link21] = lap(costMat,nonlinkMarker,0,1,noLinkCost);

                %get indices of features at time t+1 that are connected to features at time t
                indx2C = find(link21(1:m)<=n);

                %get indices of corresponding features at time t
                indx1C = link21(indx2C);

                %find the rows in "trackedFeatureNum" that are not connected to features at time t+1
                indx1U = ones(size(trackedFeatureNum,1),1);
                indx1U(indx1C) = 0;
                indx1U = find(indx1U);

                %assign space for new matrix
                tmp = zeros(size(trackedFeatureNum,1)+m-length(indx2C),t+1);

                %fill in the feature numbers at time t+1
                tmp(1:m,t+1) = [1:m]';

                %shuffle the rows from the previous times to get the correct
                %connectivity with time point t+1
                tmp(indx2C,1:t) = trackedFeatureNum(indx1C,:);

                %add rows of tracks that are not connected to points at time t+1
                tmp(max(m,length(indx1C))+1:end,1:t) = trackedFeatureNum(indx1U,:);

                %update the connectivity matrix "trackedFeatureNum"
                trackedFeatureNum = tmp;

            else %if there are no potential links

                %assign space for new matrix
                tmp = zeros(size(trackedFeatureNum,1)+m,t+1);

                %fill in the feature numbers at time t+1
                tmp(1:m,t+1) = [1:m]';

                %fill in the tracks upto time t
                tmp(m+1:end,1:t) = trackedFeatureNum;

                %update the connectivity matrix "trackedFeatureNum"
                trackedFeatureNum = tmp;

            end
                
        else %if there are no features in second frame

            %add a column of zeros for the second frame
            trackedFeatureNum = [trackedFeatureNum zeros(size(trackedFeatureNum,1),1)];

        end %(if m ~= 0 ... else ...)

    else %if there are no feature in first frame

        if m ~= 0 %if there are features in second frame

            %assign space for new matrix
            tmp = zeros(size(trackedFeatureNum,1)+m,t+1);

            %fill in the feature numbers at time t+1
            tmp(1:m,t+1) = [1:m]';

            %fill in the tracks upto time t
            tmp(m+1:end,1:t) = trackedFeatureNum;

            %update the connectivity matrix "trackedFeatureNum"
            trackedFeatureNum = tmp;

        else %if there are no features in second frame

            %add a column of zeros for the second frame
            trackedFeatureNum = [trackedFeatureNum zeros(size(trackedFeatureNum,1),1)];

        end %(if m ~= 0 ... else ...)

    end %(if n ~= 0 ... else ...)

end %(for t=1:numTimePoints-1)

%get total number of tracks
numTracks = size(trackedFeatureNum,1);

%find the time point where each track begins and then sort the vector
tpStart = zeros(numTracks,1);
for i=1:numTracks
    tpStart(i) = find((trackedFeatureNum(i,:)~=0),1,'first');
end
[tpStart,indx] = sort(tpStart);

%rearrange "trackedFeatureNum" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the 
%same time in descending order from longest to shortest.
trackedFeatureNum = trackedFeatureNum(indx,:);

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as [x y z a dx dy dz da] in image coordinate system
numRows = size(trackedFeatureNum,1);
trackedFeatureInfo = NaN*ones(numRows,8*numTimePoints);
for t=1:numTimePoints
    indx1 = find(trackedFeatureNum(:,t)~=0);
    if ~isempty(indx1) %if there are detected features in this frame
        indx2 = trackedFeatureNum(indx1,t);
        trackedFeatureInfo(indx1,8*(t-1)+1:8*t) = [movieInfo(t).xCoord(indx2,1) ...
            movieInfo(t).yCoord(indx2,1) movieInfo(t).zCoord(indx2,1) ...
            movieInfo(t).amp(indx2,1) movieInfo(t).xCoord(indx2,2) ...
            movieInfo(t).yCoord(indx2,2) movieInfo(t).zCoord(indx2,2) ...
            movieInfo(t).amp(indx2,2)];
    end
end


%%%%% ~~ the end ~~ %%%%%
