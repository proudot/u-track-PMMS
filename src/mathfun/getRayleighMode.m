function [mode,quant,grad]= getRayleighMode(data,c)
%GETRAYLEIGHMODE returns mode and an (1-c) quantile of a Rayleigh distributed data buffer
%
% SYNOPSIS [mode,quant,grad]= getRayleighMode(data,c)
%
% INPUT data : data buffer
%       c    : risk probability on which the quantile is set (optional)
%              default 0.01
%
% OUTPUT mode   : distribution mode
%        quant  : 1-c quantile
%        grad   : buffer with the gradients (allows one to display the 
%                 approximately Rayleigh distributed buffer)
%


if(nargin == 1)
   c = 0.01;
end;

if(~(isa(data,'double')))
   error('Invalid data type entered');
end;

[mode,quant,grad]=mexRayleighMode(data,c);