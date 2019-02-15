function [meanDisp,meanNNDist,fracNNDistLess2MD,fracNNDistLessMD,...
    meanPotLinksPerFeat,fracFeatsMore1PotLink,meanPotLinksPerTrack,...
    fracTracksMore1PotLink] = compDispNNDistLinks(tracks,...
    numPotLinksPerFeat,numPotLinksPerTrack)

%% Input

if nargin < 2 || isempty(numPotLinksPerFeat)
    numPotLinksPerFeat = [];
end

if nargin < 3 || isempty(numPorLinksPerTrack)
    numPotLinksPerTrack = [];
end

%% general information

%get number of tracks and number of frames
if isstruct(tracks) %if tracks are in structure format
    numTracks = length(tracks);
    tmp = vertcat(tracks.seqOfEvents);
    numFrames = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numFrames] = size(tracks);
    numFrames = numFrames/8;
end

%put tracks in matrix format if necessary
if isstruct(tracks)

    inputStructure = tracks;
    clear tracks
    
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end

    trackStartRow = ones(numTracks,1);
    for iTrack = 2 : numTracks
        trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
    end

    tracks = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numFrames);
    for i = 1 : numTracks
        startTime = inputStructure(i).seqOfEvents(1,1);
        endTime   = inputStructure(i).seqOfEvents(end,1);
        tracks(trackStartRow(i):trackStartRow(i)+...
            numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
            inputStructure(i).tracksCoordAmpCG;
    end
    
end

%% possible links per feature

if ~isempty(numPotLinksPerFeat)
    meanPotLinksPerFeat = mean(numPotLinksPerFeat);
    fracFeatsMore1PotLink = length(find(numPotLinksPerFeat>1))/length(numPotLinksPerFeat);
else
    meanPotLinksPerFeat = [];
    fracFeatsMore1PotLink = [];
end

%% possible links per track

if ~isempty(numPotLinksPerTrack)
    meanPotLinksPerTrack = mean(numPotLinksPerTrack);
    fracTracksMore1PotLink = length(find(numPotLinksPerTrack>1))/length(numPotLinksPerTrack);
else
    meanPotLinksPerTrack = [];
    fracTracksMore1PotLink = [];
end

%% mean displacement

xCoord = tracks(:,1:8:end);
yCoord = tracks(:,2:8:end);

xDisp = xCoord(:,2:end) - xCoord(:,1:end-1);
yDisp = yCoord(:,2:end) - yCoord(:,1:end-1);

xDisp = xDisp(:);
yDisp = yDisp(:);

xDisp = abs(xDisp(~isnan(xDisp)));
yDisp = abs(yDisp(~isnan(yDisp)));

meanDisp = mean(sqrt(xDisp.^2+yDisp.^2));

%% nearest-neighbor distances

nnDist = [];
for iFrame = 1 : numFrames
    
    xCoord1 = xCoord(:,iFrame);
    xCoord1 = xCoord1(~isnan(xCoord1));
    yCoord1 = yCoord(:,iFrame);
    yCoord1 = yCoord1(~isnan(yCoord1));
    
    if length(xCoord1) > 1
        
        featureDist = createDistanceMatrix([xCoord1 yCoord1],...
            [xCoord1 yCoord1]);
        
        featureDist = sort(featureDist,2);
        featureDist = featureDist(:,2);
        nnDist = [nnDist; featureDist]; %#ok<AGROW>
        
    end
    
end

meanNNDist = mean(nnDist);

fracNNDistLessMD = length(find(nnDist<meanDisp))/length(nnDist);
fracNNDistLess2MD = length(find(nnDist<(2*meanDisp)))/length(nnDist);

%% ~~~ the end ~~~
