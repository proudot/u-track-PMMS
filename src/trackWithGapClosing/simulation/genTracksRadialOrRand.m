function tracks = genTracksRadialOrRand(numTracks,radialArrangement)

%generate track starts, assuming that center is at 0,0
startAngle = rand(numTracks,1)*2*pi;
startPos = rand(numTracks,1)*5 + 5;
xCoordStart = startPos.*cos(startAngle);
yCoordStart = startPos.*sin(startAngle);

%generate track directions based on whether arrangement is radial or random
if radialArrangement
    deltaAngle = startAngle;
else
    deltaAngle = rand(numTracks,1)*2*pi;
end

%generate track ends
xCoordEnd = xCoordStart + cos(deltaAngle);
yCoordEnd = yCoordStart + sin(deltaAngle);

%put tracks together for output
tracks = repmat(struct('tracksCoordAmpCG',[],'tracksFeatIndxCG',[],...
    'seqOfEvents',[]),numTracks,1);
for iTrack = 1 : numTracks
    tracks(iTrack).tracksCoordAmpCG = ...
        [xCoordStart(iTrack) yCoordStart(iTrack) zeros(1,6) ...
        xCoordEnd(iTrack) yCoordEnd(iTrack) zeros(1,6)];
    tracks(iTrack).tracksFeatIndxCG = iTrack*ones(1,2);
    tracks(iTrack).seqOfEvents = [1 1 1 NaN; 2 2 1 NaN];
end

