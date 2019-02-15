%TRACKOBJ is a wrapper for trackCloseGapsKalman
%
% CONSTRUCTOR: obj = trackObj(coords,parameters,options)
%   IN : coords: variable containing the necessary information to construct
%                the inut argument movieInfo of trackCloseGapsKalman. Set
%                options to select the input parser.
%		parameters: (opt) Structure with fields
%                     .default : selection of default settings
%                        'standard' : standard defaults
%                        'brownian' : brownian motion only
%                     .costMatrices; gapCloseParm; kalmanFunctions
%                        optional fields to modify the default settings
%                    If parameters is not supplied, default 'brownian' is
%                      chosen
%		options: (opt) structure with fields
%                      .inputParser : selection of how to treat input
%                         'movieInfo' - data is already supplied as
%                            movieInfo. Parser will add standard deviations
%                            if necessary. This is the default selection if
%                            options is not supplied.
%                            inputInfo: struct with fields
%                             .fieldNames: 1-by-4 cell array of field names that
%                                should become 'xCoord', 'yCoord', 'zCoord', and
%                                'amp', respectively. Leave empty where you have no data
%                                (leaving zCoord empty will make the problem 2D)
%                                NOTE: use double-curly brackets to put a
%                                cell array into a field of a structure.
%                              .correctCentroid : 0: no correction
%                                (default), 1: correct translation, 2;
%                                correct translation and rotation
%                          'initCoord' - data is in initCoord-form, i.e.
%                             there is a field 'allCoord' with
%                             [x,y,z,sx,sy,sz], and a field amp
%                             inputInfo: struct with fields
%                               .amp - field name to use for amp. Default:
%                                 'amp'
%                               .correctCentroid -  0: no correction
%                                (default), 1: correct translation, 2:
%                                correct translation and rotation
%                           'connectiveComponent' - data is connective
%                                component as used for cexSegment with structure
%                                connComp(timepoints) = struct('Connectivity',
%                                   [], 'ImageSize', [], 'NumObjects', [], 
%                                   'PixelIdxList', []);
%                      .inputInfo  : additional information for the parser
%                      .outputStyle : settings for outputStyle (see
%                          property description of outputStyle)
%                      .plotOptions : settings for plot options (see
%                          property description of plotOptions)
%                      .name : name of data set. Default ''
%                      .additionalOutput : if 1, kalmanInfoLink and errFlag
%                          are also stored, if 0, they are discarded.
%                          Default: 0;
%
%   OUT: obj: trackObject with the following properties and methods:
%
% PROPERTIES
%   movieInfo/costMatrices/gapCloseParam/kalmanFunctions/probDim
%       input for trackCloseGapsKalman
%   rawMovieInfo is the movieInfo as supplied by the user (i.e. original
%       positions)
%   tracksFinal/kalmanInfoLink/errFlag
%       output of trackCloseGapsKalman
%   results : results according to the outputStyle.
%           If outputStyle is:
%             'standard': obj.results == obj.tracksFinal
%             'raw' : tracksFinal with uncorrected coordinates. Default,
%                   since 'standard' output can easily be obtained from
%                   obj.tracksFinal
%   outputStyle : structure with fields
%       .name : name of the style
%       .options : additional options
%   plotOptions : structure with fields
%       timeRange,colorTime,markerType,indicateSE,axH,image,flipXY
%           image, flipXY will only be considered for 2D data
%           timeRange    : 2-element row vector indicating time range
%                           to plot.  Optional. Default: whole movie.
%           colorTime    : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                           blue, red, etc.
%                           Optional. Default: 'k'.
%           markerType   : String indicating marker type for plotting.
%                           Only used if colorTime is not '1'.
%                           Optional. Default: 'none'.
%           indicateSE   : 1 if track starts and ends are to be indicated
%                           with circles and squares, respectively; 0
%                           otherwise. Optional. Default: 1.
%           axH          : if non-empty, data will be plotted into axes
%                           specified by axH. If 1, a new figure will be
%                           opened. If 0, plot will be into current axes.
%                           Default: 0
%           image        : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image
%           flipXY       : 1 if x and y coord should be flipped for
%                           plotting. Optional. Default: 0.
%           useRaw       : 1 if raw (unaligned) coordinates should be used
%                           for plotting. Optional. Default: 0;
%           beforeAfter  : 2-element vector with # of timepoints before and
%                           after current timepoint that is to be plotted.
%                           Used with method plotFrame
%
%    name : name of data set
%
% METHODS
%
%   trackObj = trackObj(coords,parameters,options)
%       Constructor (see above)
%
%   [tracksFinal,kalmanInfoLink,errFlag] = run(trackObj,verbose,saveResults)
%       This is the main tracking function
%       verbose      : 1 to show calculation progress, 0 otherwise.
%                      Optional. Default: 1.
%       saveResults  : 0 if no saving is requested. Default
%                      If saving is requested, structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                      Or []. Default: trackedFeatures in directory
%                      where run is initiated.
%                      Whole structure optional.
%       For description of output, please see tracksCloseGapsKalman
%
%    plot(trackObj,plotOpt)
%       Plots the track result. Uses options defined in obj.plotOptions,
%       plotOpt is a structure of the same form as plotOptions and allows
%       to overrule plotOptions.
%
%    plotFrame(trackObj,t,ah)
%       plots the track around frame t into the axes specified by ah
%
% REMARKS trackObj is a handle class. Thus, it will be passed by reference!
%         Note this is still work in progress.
%
% created with MATLAB ver.: 7.7.0.2162 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 14-Aug-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef trackObj<handle
    properties
        % input
        movieInfo
        costMatrices
        gapCloseParam
        kalmanFunctions
        probDim
        % output
        tracksFinal = [];
        kalmanInfoLink = [];
        errFlag = [];
        % additional properties
        outputStyle = struct('name','raw','options','');
        plotOptions = struct('timeRange',[],'colorTime','1',...
            'markerType','none','indicateSE',1,'axH',0,...
            'image',[],'flipXY',0,'useRaw',0,'beforeAfter',[2,1]);
        name = '';
        additionalOutput  = false;
        nTimepoints = 0;
        rawMovieInfo % raw coordinates (before alignment)
        
    end % static props
    properties (Dependent)
        results % property depending on outputStyle
    end % dependent properties
    properties (Hidden)
        % set here default options for outputStyles. For every name, there
        % need to be defaults for all corresponding options
        outputStyleDefOpt = struct('standard',[],...
            'raw',[]...
            );
        % defaults for inputInfo
        % correctCentroid: 0 - none; 1- only centroid; 2 - centroid and rot
        % findNames: check for fieldnames
        inputInfoDefOpt = struct('movieInfo',...
            struct('fieldNames',{{'xCoord','yCoord','zCoord','amp'}},...
            'correctCentroid',0,'findNames',0),...
            'initCoord',...
            struct('correctCentroid',0));
        
    end % hidden properties
    methods
        %================
        %% CONSTRUCTOR
        %================
        function obj = trackObj(coords,parameters,options)
            % allow empty object, in case we want to run a batch of tracks
            if nargin < 1 || isempty(coords)
                obj.movieInfo = [];
                return
            end
            
            % set defaults
            if nargin < 2 || isempty(parameters) || ~isfield(parameters,'default')
                parameters.default = 'brownian';
            end
            if nargin < 3 || isempty(options) || ~isfield(options,'inputParser')
                options.inputParser = 'movieInfo';
            end
            if ~isfield(options,'inputInfo')
                options.inputInfo = [];
            end
            
            % set options
            options = obj.setOptions(options);
            
            % parse input to set movieInfo. Return movieInfo explicitly so
            % that we can call parseInputs from itself, if necessary
            obj.movieInfo = obj.parseInputs(coords,options);
            
            % set parameters
            obj.setParameters(parameters);
            
            
            
        end
        %===================
        %% RUN
        %===================
        function [tracksFinal,kalmanInfoLink,errFlag] = run(obj,verbose,saveResults)
            
            %---- test input
            nObj = numel(obj);
            if nObj > 1 && nargout > 0
                error('cannot return output arguments with multiple objects')
            end
            if nargin < 2 || isempty(verbose);
                verbose = 1;
            end
            if nargin < 3 || isempty(saveResults)
                saveResults = zeros(nObj,1);
            else
                nSaveResults = length(saveResults);
                if nSaveResults == 1 && nObj>1
                    saveResults = repmat(saveResults,nObj,1);
                end
                if nSaveResults ~= nObj
                    error('for multiple objects, multiple save options must be provided')
                end
            end
            
            %---- track
            for iObj = 1:nObj
                [obj(iObj).tracksFinal,obj(iObj).kalmanInfoLink,obj(iObj).errFlag]...
                    = trackCloseGapsKalman(obj(iObj).movieInfo,...
                    obj(iObj).costMatrices,obj(iObj).gapCloseParam,...
                    obj(iObj).kalmanFunctions,obj(iObj).probDim,saveResults(iObj),verbose);
            end
            
            %---- return results if requested
            if nargout > 0
                % nObj equals to 1
                tracksFinal = obj.tracksFinal;
                kalmanInfoLink = obj.kalmanInfoLink;
                errFlag = obj.errFlag;
            end
            % discard kalmanInfoLink if unnecessary
            for iObj = 1:nObj
                if obj(iObj).additionalOutput == 0
                    obj(iObj).kalmanInfoLink = [];
                    obj(iObj).errFlag = [];
                end
            end
            
        end % run
        %======================
        %% PLOT
        %======================
        function plot(obj,plotOpt)
            % todo: set callback (don't forget to set userData in
            % plotTracks)
            
            
            % loop through objects and plot
            nObj = numel(obj);
            for iObj = 1:nObj
                
                % read plotOptions individually for each object
                if nargin < 2 || isempty(plotOpt)
                    plotOptions = obj(iObj).plotOptions; %#ok<*PROP>
                else
                    plotOptions = obj(iObj).plotOptions;
                    % overwrite plotOptions
                    for fn = fieldnames(plotOpt)'
                        plotOptions.(fn{1}) = plotOpt.(fn{1});
                    end
                end
                
                
                tmpAxH = plotOptions.axH;
                if tmpAxH == 0 || (nObj>1 && tmpAxH == 1);
                    % create figure using the name that has been provided
                    if ~isempty(obj(iObj).name)
                        figure('Name',obj(iObj).name);
                    else
                        figure('Name',sprintf('Data set %i',iObj));
                    end
                    plotOptions.plotOptions.axH = 0;
                    axes(gca);
                elseif ishandle(tmpAxH) && strcmp(get(tmpAxH,'type'),'axes')
                    plotOptions.axH = 0;
                    axes(tmpAxH);
                end
                switch obj(iObj).probDim
                    case 2
                        plotTracks2D(obj(iObj).results,...
                            plotOptions.timeRange,...
                            plotOptions.colorTime,...
                            plotOptions.markerType,...
                            plotOptions.indicateSE,...
                            plotOptions.axH,...
                            plotOptions.image,...
                            plotOptions.flipXY,0);
                    case 3
                        plotTracks3D(obj(iObj).results,...
                            plotOptions.timeRange,...
                            plotOptions.colorTime,...
                            plotOptions.markerType,...
                            plotOptions.indicateSE,...
                            plotOptions.axH);
                end
                % update plotOptions with axH
                obj(iObj).plotOptions.axH = tmpAxH;
            end
        end % plot
        %==========================
        %% PLOTFRAME
        %==========================
        function plotFrame(obj,t,axH,plotOpt,offset)
            % plotFrame plots the tracks into a GUI window
            % plotFrame(obj,t,axH,plotOpt)
            %
            % plotFrame makes the following assumptions:
            % - 2D plotting
            % - feature is getting plotted by another routine
            % - only one object
            
            if nargin < 3 || isempty(t) || isempty(axH)
                error('please supply time and axes handle for plotting the current frame')
            end
            
            % read plotOptions
            if nargin < 4 || isempty(plotOpt)
                plotOptions = obj.plotOptions;
            else
                plotOptions = obj.plotOptions;
                % overwrite plotOptions
                for fn = fieldnames(plotOpt)'
                    plotOptions.(fn{1}) = plotOpt.(fn{1});
                end
            end
            if nargin < 5 || isempty(offset)
                offset = [0,0];
            end
            
            % get start/end times
            minMaxTime = [1,length(obj.movieInfo)];
            startEndTime = NaN(1,2);
            startEndTime(1) = max(minMaxTime(1),t-plotOptions.beforeAfter(1));
            startEndTime(2) = min(minMaxTime(2),t+plotOptions.beforeAfter(2));
            
            % make sure we're getting raw coordinates
            oldStyle = obj.outputStyle.name;
            obj.outputStyle.name = 'raw';
            
            % plot
            plotTracks2D(obj.results,...
                startEndTime,...
                plotOptions.colorTime,...
                plotOptions.markerType,...
                plotOptions.indicateSE,...
                axH,...
                [],...
                plotOptions.flipXY,0,offset);
            
            % reset style
            obj.outputStyle.name = oldStyle;
            
        end
        %==========================
        %% GET.RESULTS
        %==========================
        % if there is a list of objects, return multiple outputs, as we
        % would expect of, e.g. a structure
        function varargout = get.results(obj)
            % count objects
            nObj = numel(obj);
            % assign default output
            if nargout > 1
                [varargout{1:nObj}] = deal(obj.tracksFinal);
            else
                varargout{1} = obj.tracksFinal;
            end
            
            % loop objects.
            for iObj = 1:numel(obj)
                %Since obj is passed by reference, we can copy
                % to make our lives easier
                cObj = obj(iObj);
                % only read if there is anything to do
                if ~isempty(cObj.tracksFinal)
                    % switch according to outputStyle
                    switch cObj.outputStyle.name
                        case 'standard'
                            % done already
                        case 'raw'
                            % replace coordinates with raw coordinates
                            tracks = cObj.tracksFinal;
                            
                            % read rawMovieInfo so that we do not keep
                            % calling the getMethod inside a long loop
                            rawMovieInfo = cObj.rawMovieInfo;
                            
                            % copied from makiGenerateTracks
                            for iTrack = 1 : length(tracks)
                                
                                nCompoundTracks =  size(tracks(iTrack).tracksFeatIndxCG,1);
                                
                                if nCompoundTracks == 1
                                    % single track. Everything is easy
                                    
                                    %fetch the start and end time of this track
                                    startTime = tracks(iTrack).seqOfEvents(1,1);
                                    endTime = tracks(iTrack).seqOfEvents(2,1);
                                    
                                else
                                    % tracks will start early and end late
                                    startTime = min(tracks(iTrack).seqOfEvents(:,1));
                                    endTime = max(tracks(iTrack).seqOfEvents(:,1));
                                end
                                
                                %go over all frames where this track exists
                                for iFrame =  startTime : endTime
                                    
                                    % loop through compound tracks
                                    for iCT = 1:nCompoundTracks
                                        
                                        %get the feature making up this track in this frame
                                        iFeature = tracks(iTrack).tracksFeatIndxCG(iCT,iFrame-startTime+1);
                                        
                                        %if there is a feature (not a gap)
                                        if iFeature ~= 0
                                            
                                            switch cObj.probDim
                                                case 2
                                                    %replace coordiantes and their stds
                                                    tracks(iTrack).tracksCoordAmpCG(iCT,(iFrame-startTime)*8+1:...
                                                        (iFrame-startTime)*8+2) = [rawMovieInfo(iFrame).xCoord(iFeature,1),...
                                                        rawMovieInfo(iFrame).yCoord(iFeature,1)];
                                                    tracks(iTrack).tracksCoordAmpCG(iCT,(iFrame-startTime)*8+5:...
                                                        (iFrame-startTime)*8+6) = [rawMovieInfo(iFrame).xCoord(iFeature,2),...
                                                        rawMovieInfo(iFrame).yCoord(iFeature,2)];
                                                case 3
                                                    %replace coordiantes and their stds
                                                    tracks(iTrack).tracksCoordAmpCG(iCT,(iFrame-startTime)*8+1:...
                                                        (iFrame-startTime)*8+3) = [rawMovieInfo(iFrame).xCoord(iFeature,1),...
                                                        rawMovieInfo(iFrame).yCoord(iFeature,1),...
                                                        rawMovieInfo(iFrame).zCoord(iFeature,1)];
                                                    tracks(iTrack).tracksCoordAmpCG(iCT,(iFrame-startTime)*8+5:...
                                                        (iFrame-startTime)*8+7) = [rawMovieInfo(iFrame).xCoord(iFeature,2),...
                                                        rawMovieInfo(iFrame).yCoord(iFeature,2),...
                                                        rawMovieInfo(iFrame).zCoord(iFeature,2)];
                                                    
                                                otherwise
                                                    error('dimensionality %i not implemented yet',cObj.probDim)
                                            end
                                            
                                        end
                                    end % loop compound tracks
                                    
                                end
                                
                            end %(for iTrack = 1 : numTracks)
                            
                            % assign tracks to output
                            varargout{iObj} = tracks;
                            
                        otherwise
                            error('unrecognized output style %s',cObj.outputStyle.name)
                    end
                end
            end
        end
        %============
        %% INSPECT
        %============
        function stats = inspect(obj,useRaw)
            % reads data from movieInfo to show distribution of nearest neighbor displacements etc. to allow better choice of parameters
            % stats.nnDisplacement
            % stats.nnDispAmp
            % stats.nnDistance
            
            if nargin < 2 || isempty(useRaw)
                useRaw = false;
            end
            
            stats = struct('nnDisplacement',{cell(obj.nTimepoints-1,2)},...
                'nnDispAmp',{cell(obj.nTimepoints-1,2)},...
                'nnDistance',{cell(obj.nTimepoints)},...
                'nnAmp',{cell(obj.nTimepoints)});
            
            goodTimes = arrayfun(@(x)(~isempty(x.xCoord)),obj.movieInfo);
            
            % it is possible to do this all without loops, but with lots of
            % features, there may be memory problems.
            
            %         % get nearest neighbor distance. Note that the minimum in the
            %         % distance matrix is all zeros
            %         nnd = arrayfun(@(x)(createDistanceMatrix([x.xCoord(:,1),x.yCoord(:,1),x.zCoord(:,1)],...
            %             [x.xCoord(:,1),x.yCoord(:,1),x.zCoord(:,1)])),obj.movieInfo(goodTimes),...
            %             'UniformOutput',false);
            %         ss = cellfun(@(x)(sort(x,2,'ascend')),nnd,'UniformOutput',false);
            %         stats.nnDistance(goodTimes) = cellfun(@(x)(x(:,2)),ss,'UniformOutput',false);
            %         clear nnd ss
            %         % get nearest neighbor amplitude
            %         nnd = arrayfun(@(x)(createDistanceMatrix(x.amp(:,1),x.amp(:,1))),...
            %             obj.movieInfo(goodTimes),'UniformOutput',false);
            %         ss = cellfun(@(x)(sort(x,2,'ascend')),nnd,'UniformOutput',false);
            %         stats.nnAmp(goodTimes) = cellfun(@(x)(x(:,2)),ss,'UniformOutput',false);
            %         clear nnd ss
            %
            %         goodPairs = find(goodTimes(1:end-1) & goodTimes(2:end));
            %         nnd = arrayfun(@(x,y)(createDistanceMatrix(...
            %                          [x.xCoord(:,1),x.yCoord(:,1),x.zCoord(:,1)],...
            %                          [y.xCoord(:,1),y.yCoord(:,1),y.zCoord(:,1)])),...
            %                          obj.movieInfo(goodPairs),...
            %                          obj.movieInfo(goodPairs+1),...
            %                          'UniformOutput',false);
            
            if useRaw
                movieInfo = obj.rawMovieInfo;
            else
                movieInfo = obj.movieInfo;
            end
            
            
            for t = find(goodTimes(:))'
                
                % get nn distance
                switch obj.probDim
                    case 2
                        xyz = [movieInfo(t).xCoord(:,1),movieInfo(t).yCoord(:,1)];
                        xyz(end,3) = 0;
                    case 3
                        xyz = [movieInfo(t).xCoord(:,1),movieInfo(t).yCoord(:,1),movieInfo(t).zCoord(:,1)];
                    otherwise
                        error('dimensionality %i not supported yet',obj.probDim)
                end
                if size(xyz,1) > 1
                nnd = createDistanceMatrix(xyz,xyz);
                nnd = sort(nnd,2);
                stats.nnDistance{t} = nnd(:,2); % do not get zero-diagonal
                
                % get nn amplitude
                nnd = createDistanceMatrix(movieInfo(t).amp(:,1),movieInfo(t).amp(:,1));
                nnd = sort(nnd,2);
                stats.nnAmp{t} = nnd(:,2); % do not get zero-diagonal
                end % check for at least two spots
                if t < obj.nTimepoints && all(goodTimes(t:t+1))
                    % get nnDisplacement
                    switch obj.probDim
                        case 2
                            nnd = createDistanceMatrix(xyz(:,1:2),...
                                [movieInfo(t+1).xCoord(:,1),movieInfo(t+1).yCoord(:,1)]);
                        case 3
                            nnd = createDistanceMatrix(xyz,...
                                [movieInfo(t+1).xCoord(:,1),movieInfo(t+1).yCoord(:,1),movieInfo(t+1).zCoord(:,1)]);
                        otherwise
                            error('dimensionality %i not supported yet',obj.probDim)
                    end
                    nnd = sort(nnd,2);
                    stats.nnDisplacement{t,1} = nnd(:,1);
                    if size(nnd,2) > 1
                    stats.nnDisplacement{t,2} = nnd(:,2);
                    end
                    
                    % get nnDispAmp. Careful, deltaAmp can be both positive
                    % and negative!
                    nnd = createDistanceMatrix(movieInfo(t).amp(:,1),movieInfo(t+1).amp(:,1));
                    [tmp,colIdx] = sort(abs(nnd),2);
                    rowIdx = repmat((1:size(tmp,1))',1,size(tmp,2));
                    linIdx = sub2ind(size(tmp),rowIdx(:),colIdx(:));
                    nnd(:) = nnd(linIdx);
                    stats.nnDispAmp{t,1} = nnd(:,1);
                    if size(nnd,2) > 1
                    stats.nnDispAmp{t,2} = nnd(:,2);
                    end
                end
                
            end % loop time
            
            % Plot
            figure('Name','nearest neighbour displacement')
            hc = distributionPlot(stats.nnDisplacement(:,1),[],[],[],[],[],[],1);
            figure('Name','second nearest neighbour displacement')
            hc2 = distributionPlot(stats.nnDisplacement(:,2),[],[],[],[],[],[],1);
            % read axes limits for hc2, set for hc
            y1 = ylim(hc{3});
            y2 = ylim(hc2{3});
            ylim(hc{3},[y1(1),y2(2)]);
            ylim(hc2{3},[y1(1),y2(2)]);
            figure('Name','nearest neighbour distance')
            distributionPlot(stats.nnDistance(:,1),[],[],[],[],[],[],1);
            figure('Name','nearest neighbour amplitude');
            hc = distributionPlot(stats.nnDispAmp(:,1),[],[],[],[],[],[],1);
            figure('Name','second nearest neighbour amplitude');
            hc2 = distributionPlot(stats.nnDispAmp(:,2),[],[],[],[],[],[],1);
            y1 = ylim(hc{3});
            y2 = ylim(hc2{3});
            ylim(hc{3},[y1(1),y2(2)]);
            ylim(hc2{3},[y1(1),y2(2)]);
            
        end
        %% INSPECTTRACKS
        function [out,figureHandles,axHandles] = inspectTracks(obj,moreData,verbose)
            %INSPECTTRACKS plots some hopefully useful data about tracks.
            % moreData: struct with length nTimepoints and fields with an entry for every feature
            % verbose : (opt) 0: no figures shown, stats calculated. 
            %           1: figures shown, stats calculated if requested. 
            %           2: figures shown, no stats calculated.
            %           Default: 0 if output is requested, 2 if no output
            %           is requested; 
            %
            %
            % plots: For every stat make a plot with all tracks in subplots
            
            
            % check for output
            if nargout > 0
                out = [];
                calcStats = true;
            else
                calcStats = false;
            end
            
            % check verbosity
            if nargin < 3 || isempty(verbose)
                verbose = 1-calcStats;
            end
            
            % update calcStats if necessary
            if verbose == 2
                calcStats = false;
            end
            
            
            % get tracksFinal
            tracksFinal = obj.tracksFinal;
            % get number of tracks
            nTracks = length(tracksFinal);
            
           
            
            % find moreData - names
            if nargin < 2 || isempty(moreData)
                nNames = 0;
                fn = {};
            else
                fn = fieldnames(moreData);
                nNames = length(fn);
            end
            
            %% create figures
            nBase = 2; % one figure for amp, one figure for displacement
            nFigures = nBase + nNames;
            figureHandles = zeros(nFigures,1);
            axHandles = zeros(nFigures,nTracks);
            
            nameList = [{'amp';'displacement'};fn];
            
             % init output. TrackStats is a struct of length nTracks with
            % fields
            % displacement, amp, fieldnames(moreData) 
            % nSegments, soe (=seqOfEvents)
            % offset : timepoint for timeIdx 0 (i.e. the first entry of a
            %   track (timeIdx=1) may be at timepoint other than 1)
            % mergeSplitGapIdx: nSegments-by-3 cell array.
            %   merge/split is [tms;otherSegment], gap is [tstart,tend]
            % info : nSegments-by-6 array with
            %   start, end, nMerges, nSplits, nGaps, nDataPoints
            %   nDataPoints is end-start+1-totalGapLength
            %
            % NOTE both amp and displacement calculations are filling some
            % of the statistics. 
            if calcStats
            nameListTS = [nameList,cell(size(nameList));{'nSegments',[];'soe',[];...
                'offset',[];'mergeSplitGapIdx',[];'info',[]}]';
            trackStats(1:nTracks,1) = struct(nameListTS{:});
            end
            
            % loop through the tracks and names and create subplot-axes that are hold
            % on already
            nRows = floor(sqrt(nTracks));
            nCols = ceil(nTracks/nRows);
            if verbose
            for f = 1:nFigures
                figureHandles(f) = figure('Name',sprintf('%s %s',obj.name,nameList{f}));
                for iTrack = 1:nTracks
                    axHandles(f,iTrack) = subplot(nRows,nCols,iTrack,'NextPlot','add');
                end
                allowaxestogrow(figureHandles(f));
            end
            end
            
            %% read and plot data
            
            for iTrack = 1:nTracks
                % find number of segments and fill other stats
                nSegments = size(tracksFinal(iTrack).tracksFeatIndxCG,1);
                
                soe = tracksFinal(iTrack).seqOfEvents;
                
                % offset is not necessarily soe(1), but is the earliest
                % start of any segment
                timeOffset = min(soe(:,1)) - 1;
                
                if calcStats
                trackStats(iTrack).nSegments = nSegments;
                trackStats(iTrack).soe = soe;
                trackStats(iTrack).offset = timeOffset;
                
                % find merge/split
                trackStats(iTrack).mergeSplitGapIdx = cell(nSegments,3);
                trackStats(iTrack).info = NaN(nSegments,6);
                mergeIdx = find(soe(:,2) == 2 & ~isnan(soe(:,4)));
                for mg = mergeIdx'
                    % there is a merge for segment mg, and a merge for the
                    % segment it points to. Keep KJ's definition of
                    % merge-time
                    trackStats(iTrack).mergeSplitGapIdx{soe(mg,3),1}(:,end+1) = [soe(mg,1);soe(mg,4)];
                    trackStats(iTrack).mergeSplitGapIdx{soe(mg,4),1}(:,end+1) = [soe(mg,1);soe(mg,3)];
                end
                splitIdx = find(soe(:,2) == 1 & ~isnan(soe(:,4)));
                for sp = splitIdx'
                    % there is a split for segment sp, and a split for the
                    % segment it points to. Keep KJ's definition of
                    % split-time
                    trackStats(iTrack).mergeSplitGapIdx{soe(sp,3),2}(:,end+1) = [soe(sp,1);soe(sp,4)];
                    trackStats(iTrack).mergeSplitGapIdx{soe(sp,4),2}(:,end+1) = [soe(sp,1);soe(sp,3)];
                end
                
                % fill number of m/s
                trackStats(iTrack).info(:,3) = ...
                    cellfun('size',trackStats(iTrack).mergeSplitGapIdx(:,1),2);
                trackStats(iTrack).info(:,4) = ...
                    cellfun('size',trackStats(iTrack).mergeSplitGapIdx(:,2),2);
                end
                
                for iFigure = 1:nFigures
                    switch nameList{iFigure}
                        case 'amp'
                            data = tracksFinal(iTrack).tracksCoordAmpCG(:,4:8:end);
                            
                            if verbose
                            plotCompTrackCore(axHandles(iFigure,iTrack),data,soe);
                            end
                            
                            % calculate additional statistics
                            if calcStats
                            trackStats(iTrack).amp = data;
                            % find datapoints
                            dataList = ~isnan(data);
                            trackStats(iTrack).info(:,6) = sum(dataList,2);
                            for iSeg = 1:nSegments
                                % if we have data at 2,3,7,9, gapIdx is 1 
                                % at [2,3], and diffIdx([2,3]) is [4,2]
                                % starts are 4,8 (3,7 + 1), ends 6,8 (3,7 +
                                % 4,2 - 1)
                                dataIdx = find(dataList(iSeg,:))';
                                diffIdx = diff(dataIdx);
                                gapIdx = diffIdx > 1;
                                if any(gapIdx)
                                trackStats(iTrack).mergeSplitGapIdx{iSeg,3} =...
                                    [dataIdx(gapIdx)+1;dataIdx(gapIdx)+diffIdx(gapIdx)-1];
                                end
                                
                                trackStats(iTrack).info(iSeg,1) = dataIdx(1);
                                trackStats(iTrack).info(iSeg,2) = dataIdx(end);
                            end
                            % count gaps
                            trackStats(iTrack).info(:,5) = cellfun('size',...
                                trackStats(iTrack).mergeSplitGapIdx(:,3),2);
                            end
                            
                        case 'displacement'
                            % read positions, copy start/end position for s/m and
                            % calculate displacement
                            [nSeg,n8] = size(tracksFinal(iTrack).tracksCoordAmpCG);
                            posIdx = sort([1:8:n8,2:8:n8,3:8:n8]);
                            pos = tracksFinal(iTrack).tracksCoordAmpCG(:,posIdx);
                            
                            % loop through segments and calculate displacement
                            displacement = NaN(nSeg,n8/8);
                            
                            for s = 1:nSeg
                                % check for merge/split
                                startIdx = soe(:,2)==1 & soe(:,3)==s;
                                endIdx = soe(:,2)==2 & soe(:,3)==s;
                                if ~isnan(soe(startIdx,4))
                                    % split - copy last position of other track
                                    % do this only if there is another
                                    % timepoint beforehand
                                    
                                    % splitTrackIdx: track from which to copy
                                    splitTrackIdx = soe(startIdx,4);
                                    otherStart = soe(soe(:,3)==splitTrackIdx & soe(:,2) == 1);
                                    
                                    if soe(startIdx,1) > otherStart
                                    % cpIdx: indices to copy
                                    cpIdx = (soe(startIdx,1)-1)*3 - (timeOffset)*3;
                                    cpIdx = cpIdx-2:cpIdx;
                                    %% cpIdx<0??
                                    pos(s,cpIdx) = pos(splitTrackIdx,cpIdx);
                                    end
                                end
                                if ~isnan(soe(endIdx,4))
                                    % merge - copy next position of other track
                                    % do this only if there is another
                                    % timepoint afterward
                                    
                                    % splitTrackIdx: track from which to copy
                                    splitTrackIdx = soe(endIdx,4);
                                    otherEnd = soe(soe(:,3)==splitTrackIdx & soe(:,2) == 2);
                                    if soe(endIdx,1) < otherEnd
                                    % cpIdx: indices to copy
                                    cpIdx = (soe(endIdx,1)+1)*3 - (timeOffset)*3;
                                    cpIdx = cpIdx-2:cpIdx;
                                    
                                    pos(s,cpIdx) = pos(splitTrackIdx,cpIdx);
                                    end
                                end
                                
                                % calculate displacement
                                segPos = reshape(pos(s,:),3,[])';
                                segPos = diff(segPos,[],1); %#ok<UDIM> % displacementVectors
                                
                                % write into displacements
                                displacement(s,:) = [normList(segPos);NaN]';
                                
                            end % loop segments
                            
                            if calcStats
                            trackStats(iTrack).displacement = displacement;
                            end
                            
                            % plot
                            if verbose
                            plotCompTrackCore(axHandles(iFigure,iTrack),displacement,soe);
                            end
                            
                        otherwise
                            idx = tracksFinal(iTrack).tracksFeatIndxCG;
                            data = NaN(size(idx));
                            for timeIdx = 1:size(idx,2)
                                goodIdx = idx(:,timeIdx) > 0;
                                data(goodIdx,timeIdx) = ...
                                    moreData(timeIdx + timeOffset).(nameList{iFigure})(idx(goodIdx,timeIdx));
                            end
                            
                            if calcStats
                            trackStats(iTrack).(nameList{iFigure}) = data;
                            end
                            if verbose
                            plotCompTrackCore(axHandles(iFigure,iTrack),data,soe);
                            end
                    end
                end
            end
            
            % set correct limits
            if verbose
            for f = 1:nFigures
                % xlim : global min/max
                xlimc = get(axHandles(f,:),'xlim');
                xlimAll = cat(1,xlimc{:});
                xmin = min(xlimAll(:,1));
                xmax = max(xlimAll(:,2));
                set(axHandles(f,:),'xlim',[xmin,xmax]);
                % ylim : global min/max
                ylimc = get(axHandles(f,:),'ylim');
                ylimAll = cat(1,ylimc{:});
                ymin = min(ylimAll(:,1));
                ymax = max(ylimAll(:,2));
                set(axHandles(f,:),'ylim',[ymin,ymax]);
            end
            end
            if calcStats
                out = trackStats;
            end
            
            
            %% subfcn plotCompTrackCore
            % from KJ's plotCompTrack.m
            
            function plotCompTrackCore(ah,valuesMatrix,seqOfEvents)
                % ah : axes handle
                % values matrix: nSegments-by-nPoints array of values to be plotted
                % seqOfEvents : track.seqOfEvents
                
                % samplingFreq: leftover from Khuloud's code
                samplingFreq = 1;
                
                %get first frame, last frame and number of frames
                firstFrame = seqOfEvents(1,1);
                lastFrame = max(seqOfEvents(:,1));
                
                %get number of segments making compound track
                numSegments = size(valuesMatrix,1);
                
                %plot values as dotted black lines, closing gaps
                for i = 1 : numSegments
                    indx = find(~isnan(valuesMatrix(i,:)));
                    plot(ah,(indx-1)*samplingFreq+firstFrame,valuesMatrix(i,indx),'k:');
                end
                
                %plot values in color, leaving gaps as blank (so that they appear as
                %dotted lines in the final figure)
                plot(ah,(firstFrame:samplingFreq:lastFrame)',valuesMatrix','marker','.');
                
                %find merges and splits
                indxSplit = (find(seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4))))';
                indxMerge = (find(seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))))';
                
                %go over all splits
                for iSplit = indxSplit
                    
                    %get time of splitting
                    timeSplit = seqOfEvents(iSplit,1);
                    
                    %determine location in valuesMatrix
                    splitLoc = (timeSplit - firstFrame) / samplingFreq + 1;
                    
                    %determine index of starting track
                    rowS = seqOfEvents(iSplit,3);
                    
                    %determine index of splitting track
                    rowSp = seqOfEvents(iSplit,4);
                    
                    %plot split as a black dash-dotted line
                    plot(ah,[timeSplit-samplingFreq timeSplit],[valuesMatrix(rowSp,splitLoc-1) ...
                        valuesMatrix(rowS,splitLoc)],'k-.')
                    
                end
                
                %go over all merges
                for iMerge = indxMerge
                    
                    %get time of merging
                    timeMerge = seqOfEvents(iMerge,1);
                    
                    %determine location in valuesMatrix
                    mergeLoc = (timeMerge - firstFrame) / samplingFreq + 1;
                    
                    %determine index of ending track
                    rowE = seqOfEvents(iMerge,3);
                    
                    %determine index of merging track
                    rowM = seqOfEvents(iMerge,4);
                    
                    %plot merge as a black dashed line
                    plot(ah,[timeMerge-samplingFreq timeMerge],[valuesMatrix(rowE,mergeLoc-1) ...
                        valuesMatrix(rowM,mergeLoc)],'k--')
                    
                end
            end
        end
        %% Get Method rawMovieInfo
        function out = get.rawMovieInfo(obj)
            % return movieInfo unless correctCoords has written
            % rawMovieInfo
            if isempty(obj.rawMovieInfo)
                out = obj.movieInfo;
            else
                out = obj.rawMovieInfo;
            end
        end
        %% setMethod gapCloseParam
        function set.gapCloseParam(obj,val)
            fn = fieldnames(val);
            for field = fn'
                field = field{1};
                
                % check for whether timeWindow exists and whether it has
                % changed
                if strcmp(field,'timeWindow') && isfield(obj.gapCloseParam,'timeWindow') && val.timeWindow ~= obj.gapCloseParam.timeWindow
                    % make sure that dependent properties are up-to-date
                    if ~isempty(obj.costMatrices) %#ok<*MCSUP> %
                        % check for the existence of fields before updating
                        if isfield(obj.costMatrices(1),'parameters') > 0 && isfield(obj.costMatrices(1).parameters,'nnWindow') %#ok<MCSUP>
                            obj.costMatrices(1).parameters.nnWindow = val.timeWindow;
                        end
                        if length(obj.costMatrices) > 1
                            if isfield(obj.costMatrices(2).parameters,'nnWindow')
                                obj.costMatrices(2).parameters.nnWindow = val.timeWindow;
                            end
                            if isfield(obj.costMatrices(2).parameters,'brownStdMult')
                                if all(obj.costMatrices(2).parameters.brownStdMult == obj.costMatrices(2).parameters.brownStdMult(1))
                                    obj.costMatrices(2).parameters.brownStdMult = repmat(obj.costMatrices(2).parameters.brownStdMult(1),val.timeWindow,1);
                                else
                                    warning('please update trackObj.costMatrices(2).parameters.brownStdMult')
                                end
                            end
                            if isfield(obj.costMatrices(2).parameters,'timeReachConfB')
                                % only update if it was the same as timeWindow before
                                if obj.costMatrices(2).parameters.timeReachConfB == obj.gapCloseParam.timeWindow
                                    obj.costMatrices(2).parameters.timeReachConfB = val.timeWindow;
                                end
                            end
                            if isfield(obj.costMatrices(2).parameters,'timeReachConfL')
                                % only update if it was the same as timeWindow before
                                if obj.costMatrices(2).parameters.timeReachConfL == obj.gapCloseParam.timeWindow
                                    obj.costMatrices(2).parameters.timeReachConfL = val.timeWindow;
                                end
                            end
                            if isfield(obj.costMatrices(2).parameters,'linStdMult')
                                if all(obj.costMatrices(2).parameters.linStdMult == obj.costMatrices(2).parameters.linStdMult(1))
                                    obj.costMatrices(2).parameters.linStdMult = repmat(obj.costMatrices(2).parameters.linStdMult(1),val.timeWindow,1);
                                else
                                    warning('please update trackObj.costMatrices(2).parameters.linStdMult')
                                end
                            end
                        end % exist costMatrices(2)
                    end % isempty costmatrices
                end % check for field timeWindow
                
                
                
                % update parameters only now so that we can compare to the
                % old timeWindow if necessary
                obj.gapCloseParam.(field) = val.(field);
            end % loop fields
            
        end % function setGapCloseParam
    end % public methods
    methods (Hidden)
        %===========================
        %% PARSE INPUTS
        %===========================
        function movieInfo = parseInputs(obj,coord,options)
            switch options.inputParser
                case 'movieInfo'
                    % set default field names
                    movieInfoFields = obj.inputInfoDefOpt.movieInfo.fieldNames;
                    if ~options.inputInfo.findNames
                        fieldList = options.inputInfo.fieldNames;
                    else
                        % if no explicit list: find fieldnames from coord
                        fieldList = fieldnames(coord);
                        nameIdx = ismember(movieInfoFields',fieldList);
                        if nameIdx(1)==0
                            error('you need at least to supply x coords')
                        end
                        fieldList = cell(size(movieInfoFields));
                        fieldList(nameIdx) = movieInfoFields(nameIdx);
                    end
                    if length(fieldList) == 3 || isempty(fieldList{3})
                        % 2D problem. remove zCoord
                        obj.probDim = 2;
                        movieInfoFields(3) = [];
                        if isempty(fieldList{3})
                            fieldList(3) = [];
                        end
                    else
                        obj.probDim = 3;
                    end
                    % create movieInfo
                    nTimepoints = length(coord);
                    tmp = movieInfoFields;
                    tmp{2,end} = [];
                    movieInfo(1:nTimepoints) = struct(tmp{:});
                    
                    % loop to fill movieInfo
                    for iName = 1:1+obj.probDim
                        for t = 1:nTimepoints
                            if ~isempty(fieldList(iName))
                                movieInfo(t).(movieInfoFields{iName}) = ...
                                    coord(t).(fieldList{iName});
                                if size(movieInfo(t).(movieInfoFields{iName}),2) == 1
                                    movieInfo(t).(movieInfoFields{iName})(:,2) = 1;
                                end
                            else
                                % set the un-supplied data to 1. This makes
                                % only really sense for amp, though
                                movieInfo(t).(movieInfoFields{iName}) = ...
                                    ones(size(movieInfo(t).xCoord));
                            end
                        end
                    end
                    
                    % correct centroid
                    correctCentroid;
                    
                case 'initCoord'
                    % check options
                    if ~isempty(options.inputInfo) && isfield(options.inputInfo,'amp')
                        ampName = options.inputInfo.amp;
                    else
                        ampName = 'amp';
                    end
                    
                    % problem dimension is automatically 3
                    obj.probDim = 3;
                    % create movieInfo
                    nTimepoints = length(coord);
                    movieInfo(1:nTimepoints) = struct('xCoord',[],...
                        'yCoord',[],'zCoord',[],'amp',[]);
                    
                    goodTimes = find(arrayfun(@(x)(~isempty(x.allCoord)),...
                        coord));
                    
                    if ~isempty(goodTimes)
                        % only do the rest if nonempty input
                        
                        for t = goodTimes'
                            movieInfo(t).xCoord = coord(t).allCoord(:,[1 4]);
                            movieInfo(t).yCoord = coord(t).allCoord(:,[2 5]);
                            movieInfo(t).zCoord = coord(t).allCoord(:,[3 6]);
                            movieInfo(t).amp = coord(t).(ampName);
                        end
                        
                        correctCentroid;
                        
                    end
   
                    
                case 'connectiveComponent'
                    % check options
                    if ~isempty(options.inputInfo) && isfield(options.inputInfo,'amp')
                        ampName = options.inputInfo.amp;
                    else
                        ampName = 'volume';
                    end
                    
                    % coords will be a table like so
                    %coords(x) = struct('Connectivity', [], ...
                    %  'ImageSize', [], 'NumObjects', [], 'PixelIdxList', []);
                    
                    voxelSize = obj.voxelSize(1) * obj.voxelSize(2) * obj.voxelSize(3);

                    
                    % problem dimension is automatically 3
                    obj.probDim = 3;
                    % create movieInfo
                    nTimepoints = length(coord);
                    movieInfo(1:nTimepoints) = struct('xCoord',[],...
                        'yCoord',[],'zCoord',[],'amp',[]);
                    
                    %goodTimes = find(arrayfun(@(x)(~isempty(x.allCoord)), coord));
                    
                    % get 'true' for all connective components that aren't empty
                    goodTimes = arrayfun(@(x)(~isempty(x)), coords);
                    % get their indices
                    %goodIndex = find(goodIndexBoolean);
                    
                    %hack
                    params.scaleSpaceFactors = [4 4 2; 8 8 4];
                    
                    if ~isempty(goodTimes)
                        % only do the rest if nonempty input
                        
                        
                        for t = goodTimes'
                           curFrame = obj.imageData.scalespace({t},params.scaleSpaceFactors);   
                           
                           for i = 1:coords.NumObjects %coords for each potential object in CC
                               % voxels in each CC object
                               [subs(:,1) subs(:,2) subs(:,3)] = ind2sub(size(curFrame), coords(t).PixelIdxList{i});
                               movieInfo(t).xCoord(i) = subs(:,1);
                               movieInfo(t).yCoord(i) = subs(:,2);
                               movieInfo(t).zCoord(i) = subs(:,3);
                               % volume
                               movieInfo(t).amp(i) = length(CC.PixelIdxList{iChrom}) * voxelSize;
                               
                               clear subs;
                           end
                        end
                        
                        correctCentroid;
                        
                    end
                    
            end % switch inputParser
            
            obj.nTimepoints = length(movieInfo);
            
            %--- nested function correctCentroid
            function correctCentroid
                % correct Centroid
                if options.inputInfo.correctCentroid > 0
                    % store rawMovieInfo
                    obj.rawMovieInfo = movieInfo;
                    
                    % if not correct rotation, co is -1, otherwise inf
                    if options.inputInfo.correctCentroid == 1;
                        co = -1;
                    else
                        co = Inf;
                    end
                    [movieInfo,obj.outputStyle.centroids,...
                        obj.outputStyle.rotMats] = alignCoords(...
                        movieInfo,{'xCoord','yCoord','zCoord','amp'},...
                        co,obj.probDim);
                end
            end
        end % parseInputs
        %===========================
        %% SET PARAMETERS
        %===========================
        function setParameters(obj,parameters)
            % read defaults
            setDefaultParameters(obj,parameters.default);
            
            % go through parameters and overwrite defaults. Blindly step
            % through fieldnames and overwrite whatever was in the input.
            % This requires that the user made sensible choices.
            for pn = fieldnames(parameters)'
                if ~strcmp(pn{1},'default')
                    for k = 1:length(parameters.(pn{1}))
                        for ppn = fieldnames(parameters.(pn{1})(k))
                            if isstruct(parameters.(pn{1})(k).(ppn{1}))
                                for pppn = fieldnames(parameters.(pn{1})(k).(ppn{1}))'
                                    try
                                        obj.(pn{1})(k).(ppn{1}).(pppn{1}) = ...
                                            parameters.(pn{1})(k).(ppn{1}).(pppn{1});
                                    catch
                                        error(...
                                            'parameters.%s(%i).%s.%s is not a recognized option',...
                                            pn{1},k,ppn{1},pppn{1})
                                    end
                                end
                            else
                                try
                                    obj.(pn{1})(k).(ppn{1}) = ...
                                        parameters.(pn{1})(k).(ppn{1});
                                catch
                                    error(...
                                        'parameters.%s(%i).%s is not a recognized option',...
                                        pn{1},k,ppn{1})
                                end
                            end
                        end
                    end
                end
            end
            
            
        end % set parameters
        %===========================
        %% SET OPTIONS
        %===========================
        function options = setOptions(obj,options)
            %		options: (opt) structure with fields
            %                      .inputParser : selection of how to treat input
            %                         'movieInfo' - data is already supplied as
            %                            movieInfo. Parser will add standard deviations
            %                            if necessary. This is the default selection if
            %                            options is not supplied.
            %                      .inputInfo  : additional information for the parser
            %                      .outputStyle : settings for outputStyle (see
            %                          property description of outputStyle)
            %                      .plotOptions : settings for plot options (see
            %                          property description of plotOptions)
            %                      .name : name of data set. Default ''
            %                      .additionalOutput : if 1, kalmanInfoLink and errFlag
            %                          are also stored, if 0, they are discarded.
            %                          Default: 0;
            % options will never be empty (there is at least the default
            % input parser. Therefore, just test for fieldnames
            
            % outputStyle
            if isfield(options,'outputStyle')
                % set defaults
                if isfield(options.outputStyle,'name')
                    obj.outputStyle.name = options.outputStyle.name;
                end
                obj.outputStyle.options = obj.outputStyleDefOpt.(...
                    obj.outputStyle.name);
                % overwrite default with additional options
                for fn = fieldnames(options.outputStyle)'
                    if strcmp(fn{1},'options')
                        for ofn = fieldnames(options.outputStyle.options)'
                            obj.outputStyle.options.(ofn{1}) = ...
                                options.outputStyle.options.(ofn{1});
                        end
                    else
                        obj.outputStyle.(fn{1}) = options.outputStyle.(fn{1});
                    end
                end
            end
            
            % inputInfo
            if isempty(options.inputInfo)
                options.inputInfo = obj.inputInfoDefOpt.(options.inputParser);
            else
                tmp = options.inputInfo;
                options.inputInfo = obj.inputInfoDefOpt.(options.inputParser);
                for fn = fieldnames(tmp)'
                    options.inputInfo.(fn{1}) = tmp.(fn{1});
                end
                if strcmp(options.inputParser,'movieInfo') && ...
                        (~isfield(tmp,'fieldNames') || isempty(tmp.fieldNames))
                    options.inputInfo.findNames = true;
                end
            end
            
            if isfield(options,'plotOptions')
                for fn = fieldnames(options.plotOptions)'
                    obj.plotOptions.(fn{1}) = options.plotOptions.(fn{1});
                end
            end
            if isfield(options,'name')
                obj.name = options.name;
            end
            if isfield(options,'additionalOutput')
                obj.additionalOutput = options.additionalOutput;
            end
        end % set options
        %===========================
        %% DEFAULTPARAMETERS
        %===========================
        function setDefaultParameters(obj,default)
            switch lower(default)
                case 'brownian'
                    % general gap closing parameters
                    obj.gapCloseParam.timeWindow = 5; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
                    obj.gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 0 otherwise.
                    obj.gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
                    
                    % cost matrix for frame-to-frame linking
                    
                    %function name
                    obj.costMatrices(1).funcName = 'costMatLinearMotionLink';
                    
                    %parameters
                    
                    parameters.linearMotion = 0; %no linear motion
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
                    parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
                    
                    parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
                    
                    obj.costMatrices(1).parameters = parameters;
                    clear parameters
                    
                    % cost matrix for gap closing
                    
                    %function name
                    obj.costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';
                    
                    %parameters
                    
                    %needed all the time
                    parameters.linearMotion = 0;%1; %use linear motion Kalman filter.
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius.
                    parameters.brownStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
                    parameters.timeReachConfB = 2; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
                    
                    parameters.ampRatioLimit = [0.5 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
                    
                    parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
                    
                    parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
                    
                    parameters.linStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
                    parameters.timeReachConfL = obj.gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
                    parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
                    
                    obj.costMatrices(2).parameters = parameters;
                    clear parameters
                    
                    % Kalman filter function names
                    
                    obj.kalmanFunctions.reserveMem = 'kalmanResMemLM';
                    obj.kalmanFunctions.initialize = 'kalmanInitLinearMotion';
                    obj.kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
                    obj.kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
                case 'standard'
                    % general gap closing parameters
                    obj.gapCloseParam.timeWindow = 5; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
                    obj.gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 0 otherwise.
                    obj.gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
                    
                    % cost matrix for frame-to-frame linking
                    
                    %function name
                    obj.costMatrices(1).funcName = 'costMatLinearMotionLink';
                    
                    %parameters
                    
                    parameters.linearMotion = 1; %use linear motion Kalman filter.
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
                    parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
                    
                    parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
                    
                    obj.costMatrices(1).parameters = parameters;
                    clear parameters
                    
                    % cost matrix for gap closing
                    
                    %function name
                    obj.costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';
                    
                    %parameters
                    
                    %needed all the time
                    parameters.linearMotion = 1;%1; %use linear motion Kalman filter.
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius.
                    parameters.brownStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
                    parameters.timeReachConfB = 2; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
                    
                    parameters.ampRatioLimit = [0.5 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
                    
                    parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
                    
                    parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
                    
                    parameters.linStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
                    parameters.timeReachConfL = obj.gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
                    parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
                    
                    obj.costMatrices(2).parameters = parameters;
                    clear parameters
                    
                    % Kalman filter function names
                    
                    obj.kalmanFunctions.reserveMem = 'kalmanResMemLM';
                    obj.kalmanFunctions.initialize = 'kalmanInitLinearMotion';
                    obj.kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
                    obj.kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
                otherwise
                    error('%s not implemented yet',default)
            end
        end
        
    end % hidden methods
    
end