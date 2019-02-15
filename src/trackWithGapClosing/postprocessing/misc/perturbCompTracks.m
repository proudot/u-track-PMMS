function tracksPerturb = perturbCompTracks(tracksFinal,fracMissing,fracMSFalse)
%PERTURBCOMPTRACKS perturbs compound tracks to mimic error in experimental data
%
%SYNOPSIS tracksPerturb = perturbCompTracks(tracksFinal,percentMissing,percentMSWrong)
%
%INPUT  tracksFinal  : Output of trackCloseGapsKalman.
%       fracMissing  : Fraction of missing particles.
%                      Optional. Default: 0.05.
%       fracMSFalse  : Fraction of merging and splitting false positives
%                      and false negatives.
%                      Optional. Default: 0.05.
%
%OUTPUT tracksPerturb: Similar to tracksFinal, but after track perturbation
%
%REMARKS Pertubration probably on the simplistic side for now.
%
%Khuloud Jaqaman, September 2013

%% Output

tracksPerturb = tracksFinal;

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--perturbCompTracks: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(fracMissing)
    fracMissing = 0.05;
end

if nargin < 3 || isempty(fracMSFalse)
    fracMSFalse = 0.05;
end

%% track perturbation

%number of compound tracks
numTracks = length(tracksFinal);

%number of frames
seqOfEvents = vertcat(tracksFinal.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%determine number of segments per compound track
numSeg = ones(numTracks,1);
for iTrack = 1 : numTracks
    numSeg(iTrack) = size(tracksFinal(iTrack).tracksFeatIndxCG,1);
end

%get the track start and end times
trackSEL = getTrackSEL(tracksFinal);

%remove fracMissing of particles, only from tracks without merging and splitting events
for iFrame = 1 : numFrames
    indxNoMSExist = find(numSeg==1&trackSEL(:,1)<iFrame&trackSEL(:,2)>iFrame);
    missIndx = randsample(indxNoMSExist,round(fracMissing*length(indxNoMSExist)));
    if ~isempty(missIndx)
        for iTrack = missIndx'
            seqOfEvents = tracksPerturb(iTrack).seqOfEvents;
            iFrameRel = iFrame - seqOfEvents(1,1) + 1;
            tracksPerturb(iTrack).tracksFeatIndxCG(iFrameRel) = 0;
            tracksPerturb(iTrack).tracksCoordAmpCG((iFrameRel-1)*8+1:iFrameRel*8) = NaN;
        end
    end
end

%remove fracMSFalse of the merging and splitting events
%for the sake of computational simplicity, replace those with false
%positive merges and splits (this keeps number of tracks the same)
indxMSAll = find(numSeg>1);
indxMS23 = find(numSeg==2|numSeg==3);
indxMSPerturb = randsample(indxMS23,round(fracMSFalse*length(indxMSAll)));
for iTrack = indxMSPerturb'
    msTime = randsample(numFrames-2,1)+1;
    eventType = randi(2);
    switch eventType
        case 1
            seqOfEvents = [[1 msTime numFrames numFrames]' [1 1 2 2]' [1 2 2 1]' NaN(4,1)];
            seqOfEvents(2,4) = 1;
            tracksFeatIndx = zeros(2,numFrames);
            tracksFeatIndx(1,:) = 1;
            tracksFeatIndx(2,msTime:end) = 1;
            tracksCoordAmp = NaN(2,numFrames*8);
            tracksCoordAmp(1,:) = 1;
            tracksCoordAmp(2,(msTime-1)*8+1:end) = 1;
        case 2
            seqOfEvents = [[1 1 msTime numFrames]' [1 1 2 2]' [1 2 2 1]' NaN(4,1)];
            seqOfEvents(3,4) = 1;
            tracksFeatIndx = zeros(2,numFrames);
            tracksFeatIndx(1,:) = 1;
            tracksFeatIndx(2,1:msTime-1) = 1;
            tracksCoordAmp = NaN(2,numFrames*8);
            tracksCoordAmp(1,:) = 1;
            tracksCoordAmp(2,1:(msTime-1)*8) = 1;
    end
    tracksPerturb(iTrack).tracksFeatIndxCG = tracksFeatIndx;
    tracksPerturb(iTrack).tracksCoordAmpCG = tracksCoordAmp;
    tracksPerturb(iTrack).seqOfEvents = seqOfEvents;
end

%%%%% ~~ the end ~~ %%%%%

