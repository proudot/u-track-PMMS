% MISC
%
% Files
%   analyzeQDotBlinking        - define imaging conditions to be compared
%   batchSfaAllEnd             - 
%   chooseTracks               - outputs the indices of tracks that satisfy the input criteria
%   compDispNNDistLinks        - Input
%   compLocalDisp2NNDist       - compares for every particle its displacement to its nearest neighbor distance
%   convertMat2Dat             - write the output of trackCloseGapsKalman into text files for the Kusumi lab
%   convertMat2Struct2         - get number of tracks
%   convJCCellsToMat           - 
%   convMetamorphTracks2Struct - get number of tracks
%   convStruct2MatIgnoreMS     - converts tracks from structure format to matrix format, ignoring merges/splits.
%   convStruct2MatNoMS         - converts tracks from structure format to matrix format, provided there are NO merges/splits.
%   convStruct2SparseMatNoMS   - converts tracks from structure format to sparse matrix format, provided there are NO merges/splits.
%   convTrackFormatDefault2Alt - converts compound tracks from default format to alternative format
%   evalGapsMergesSplits       - calculates displacements over gaps, merges and splits and compares them to linking displacements
%   findLifetimesStatusSimple  - find lifetimes and status of all tracks in the trackInfo matrix
%   findTrackGaps              - finds the gaps in each track and gives back their location and length
%   getAppearanceInfo          - extracts the appearance position of particles in a time-lapse sequence
%   getNumSegmentsPerTrack     - GETNUMSEGMENTS gives back the number of segments in each compound track
%   getTrackSEL                - outputs track start times, end times and lifetimes
%   intModalAnalysisASE        - looks for intensity distribution modes at track starts and ends and throughout track lifetimes
%   nnDistFromDetection        - 
%   perturbCompTracks          - perturbs compound tracks to mimic error in experimental data
%   scriptCollectTracks        - define batch job locations
