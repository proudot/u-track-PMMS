%FITSEGMENT2D Fit a 2-D Diffraction-limited Segment function to data in an image window.
%    [prmVect prmStd C res J] = fitSegment2D(data, prmVect, mode)
%
%    Symbols: xp : x-position
%             yp : y-position
%              A : amplitude
%              l : length of the segment
%          sigma : half width of the segment
%              t : angle
%              c : background
%
%    The origin is defined at the center of the input window.
%
%    Inputs:     data : 2-D image array
%             prmVect : parameter vector with order: [xp, yp, A, l, s, t, c]
%                mode : string that defines parameters to be optimized; any among 'xyalstc'
%
%    Outputs: prmVect : parameter vector
%              prmStd : parameter standard deviations
%                   C : covariance matrix
%                 res : residuals
%                   J : Jacobian
%
% Axis conventions: image processing, see meshgrid
%
% Example: [prmVect prmStd C res J] = fitSegment2D(data, [0 0 max(data(:)) 10 1.5 pi/6 min(data(:))], 'xyalstc');

% (c) Sylvain Berlemont, 2011