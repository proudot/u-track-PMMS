function d = dirList(dpath)
% function d = dirList(dpath)
%
% dirList is the same as 'dir' but returns only visible directories.
%
% Francois Aguet, 11/02/2009


d = dir(dpath);

% remove entries that are not directories
d([d.isdir]==0) = [];

% remove invisible directories
d(cellfun(@isInvisible, {d(:).name})) = [];

function v = isInvisible(directory)
v = strcmp(directory(1), '.');