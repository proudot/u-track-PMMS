function copyright = getLCCBCopyright()
%
% This is a user-defined function used in LCCB software. 
% It is called when any GUI is generated. It configures the copyright
% information.
%
% Input: 
%
%
% Output:
%
%   copyright - String: copyright and version information
%
% Chuangang Ren, 11/2010
% Sebastien Besson, Feb 2013

% Set year and version information
str_year = '2014';
copyright = sprintf('Copyright %s LCCB', str_year);
