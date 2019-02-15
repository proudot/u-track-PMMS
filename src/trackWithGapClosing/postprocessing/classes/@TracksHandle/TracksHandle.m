classdef TracksHandle < Tracks & dynamicprops
% TracksHandle is a Tracks implementation optimized for serving the new
% properties such as X, Y, Z while also providing backwards-compatability
% with the tracksFinal struct
%
% See also Tracks
%
% Mark Kittisopikul, January 2015
    properties
        % x coordinates as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,1:8:end)
        % Column is relative to the startFrame
        x
        % Y coordinates as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,2:8:end)
        % Column is relative to the startFrame
        y
        % Z coordinates as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,3:8:end)
        % Column is relative to the startFrame
        z
        % Amplitude as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,4:8:end)
        % Column is relative to the startFrame
        A
        % X uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,5:8:end)
        % Column is relative to the startFrame
        dx
        % Y uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,6:8:end)
        % Column is relative to the startFrame
        dy
        % Z uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,7:8:end)
        % Column is relative to the startFrame
        dz
        % A uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,8:8:end)
        % Column is relative to the startFrame
        dA
        % Column vector for the absolute time point each segment started
        segmentStartFrame
        % Column vector for the absolute time point each segmented ended
        segmentEndFrame
        % Column vector for the segment from which the segment originated
        parentSegment
        % Column vector for the segment into which the segment merged
        spouseSegment
        % Scalar absolute frame at which the compound track starts
        startFrame
        % Scalar absolute at which the compound track ends
        endFrame
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
    end
    properties (Dependent = true)
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % 2D matrix, corresponds to tracksCoordAmpCG3D(:,:)
        tracksCoordAmpCG
        % numEvents x 4 matrix. See class description.
        seqOfEvents
        % Number of segments in each compound track
        % see also getNumSegments
        numSegments
        % Number of frames in which each compound track exists
        numFrames
    end
    properties (Access = protected, Transient)
        cache
    end
    methods
        function obj = TracksHandle(tracks,movieInfo)
            % Takes a tracksFinal structure from trackCloseGapsKalman
            if(nargin ~= 0)
                if(~isstruct(tracks))
                    tracks = convertMat2Struct2(tracks);
                end
                if(nargin > 1)
                    tracks = normalizeTracks(tracks,movieInfo);
                end
                obj(numel(tracks)) = TracksHandle();
                obj = reshape(obj,size(tracks));
                [obj.tracksFeatIndxCG] = deal(tracks.tracksFeatIndxCG);
                [obj.seqOfEvents] = deal(tracks.seqOfEvents);
                [obj.tracksCoordAmpCG] = deal(tracks.tracksCoordAmpCG);
                obj.reindex();
            end
        end
        function set.tracksFeatIndxCG(obj,tracksFeatIndxCG)
            obj.tracksFeatIndxCG = tracksFeatIndxCG;
        end
        function set.tracksCoordAmpCG(obj,tracksCoordAmpCG)
            threeD = reshape(tracksCoordAmpCG,size(tracksCoordAmpCG,1),8,[]);
            obj.tracksCoordAmpCG3D = threeD;
        end
        function tracksCoordAmpCG = get.tracksCoordAmpCG(obj)
            tracksCoordAmpCG = obj.tracksCoordAmpCG3D(:,:);
        end
        function set.tracksCoordAmpCG3D(obj,tracksCoordAmpCG3D)
            obj.x  = tracksCoordAmpCG3D(:,1,:);
            obj.y  = tracksCoordAmpCG3D(:,2,:);
            obj.z  = tracksCoordAmpCG3D(:,3,:);
            obj.A  = tracksCoordAmpCG3D(:,4,:);
            obj.dx = tracksCoordAmpCG3D(:,5,:);
            obj.dy = tracksCoordAmpCG3D(:,6,:);
            obj.dz = tracksCoordAmpCG3D(:,7,:);
            obj.dA = tracksCoordAmpCG3D(:,8,:);
            
            obj.x  = obj.x(:,:);
            obj.y  = obj.y(:,:);
            obj.z  = obj.z(:,:);
            obj.A  = obj.A(:,:);
            obj.dx = obj.dx(:,:);
            obj.dy = obj.dy(:,:);
            obj.dz = obj.dz(:,:);
            obj.dA = obj.dA(:,:);
        end
        function tracksCoordAmpCG3D = get.tracksCoordAmpCG3D(obj)
            tracksCoordAmpCG3D = zeros(obj.numSegments,8,obj.numFrames);
            tracksCoordAmpCG3D(:,1,:) = obj.x;
            tracksCoordAmpCG3D(:,2,:) = obj.y;
            tracksCoordAmpCG3D(:,3,:) = obj.z;
            tracksCoordAmpCG3D(:,4,:) = obj.A;
            tracksCoordAmpCG3D(:,5,:) = obj.dx;
            tracksCoordAmpCG3D(:,6,:) = obj.dy;
            tracksCoordAmpCG3D(:,7,:) = obj.dz;
            tracksCoordAmpCG3D(:,8,:) = obj.dA;
        end
        function set.seqOfEvents(obj,seqOfEvents)
            init = zeros(obj.numSegments,1);
            
            obj.segmentStartFrame = init;
            startIdx = seqOfEvents(:,2) == 1;
            obj.segmentStartFrame(seqOfEvents(startIdx,3)) = seqOfEvents(startIdx,1);
            
            obj.segmentEndFrame = init;
            endIdx = seqOfEvents(:,2) == 2;
            obj.segmentEndFrame(seqOfEvents(endIdx,3)) = seqOfEvents(endIdx,1);
            endIdx = endIdx & ~isnan(seqOfEvents(:,4));
            obj.segmentEndFrame(seqOfEvents(endIdx,3)) = seqOfEvents(endIdx,1) - 1;
            
            obj.parentSegment = init;
            idx = seqOfEvents(:,2) == 1;
            obj.parentSegment(seqOfEvents(idx,3)) = seqOfEvents(idx,4);
            
            obj.spouseSegment = init;
            idx = seqOfEvents(:,2) == 2;
            obj.spouseSegment(seqOfEvents(idx,3)) = seqOfEvents(idx,4);
            
            obj.startFrame = seqOfEvents(1,1);
            obj.endFrame = seqOfEvents(end,1);
        end
        function seqOfEvents = get.seqOfEvents(obj)
            seqOfEvents = zeros(obj.numSegments*2,4);
            
            startIdx = 1:obj.numSegments;
            seqOfEvents(startIdx,1) = obj.segmentStartFrame;
            seqOfEvents(startIdx,2) = 1;
            seqOfEvents(startIdx,3) = 1:obj.numSegments;
            seqOfEvents(startIdx,4) = obj.parentSegment;
            
            % For segments that merge, their true endFrames are offset by 1
            E = obj.segmentEndFrame;
            segmentMerged = ~isnan(obj.spouseSegment);
            E(segmentMerged) = E(segmentMerged) + 1;
            
            endIdx = (1:obj.numSegments) + obj.numSegments;
            seqOfEvents(endIdx,1) = E;
            seqOfEvents(endIdx,2) = 2;
            seqOfEvents(endIdx,3) = 1:obj.numSegments;
            seqOfEvents(endIdx,4) = obj.spouseSegment;
            
            seqOfEvents = sortrows(seqOfEvents);
        end
        function resetCache(obj)
            c = struct();
            c.tracksCoordAmpCG = [];
            obj.cache = c;
        end
        function N = get.numSegments(obj)
            N = size(obj.x,1);
        end
        function N = get.numFrames(obj)
            N = size(obj.x,2);
        end
        function P = addprop(obj,propName)
            if(~isscalar(obj))
                P = arrayfun(@(x) addprop(x,propName),obj,'UniformOutput',false);
                P = [P{:}];
                P = reshape(P,size(obj));
            else
                P = addprop@dynamicprops(obj,propName);
            end
        end
    end
end
