
%for the initial simple linking from frame to frame
costMatrices(1).costMatFun = 'costMatSimple';
costMatrices(1).costMatParam = struct('searchRadius',3,'maxAmpRatio',2,'noLnkPrctl',-1);

%for linking between frames again using statistical data on the tracks
costMatrices(2).costMatFun = 'costMatLogL';
costMatrices(2).costMatParam = struct('cutoffProbD',0.99,...
    'cutoffProbA',0.99,'noLnkPrctl',-1);

%for gap closing
costMatrices(3).costMatFun = 'costMatCloseGaps';
costMatrices(3).costMatParam = struct('cutoffProbD1',0.99,'cutoffProbA1',0.99,...
    'cutoffProbD2',0.99,'cutoffProbA2',0.99,'noLnkPrctl',-1);

%for merging and splitting
costMatrices(4).costMatFun = 'costVecLinkMS';
costMatrices(4).costMatParam = struct('cutoffProbD',0.99,'cutoffProbA',0.99);

%gap closing parameters
gapCloseParam.timeWindow = 5;
gapCloseParam.mergeSplit = 0;
% gapCloseParam.segmentLength = 100;

%iteration parameters
iterParam.tolerance = 0.05;
iterParam.lenFrac = 0.5;

