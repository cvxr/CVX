function s = cvx_where

%CVX_WHERE    Returns the location of the CVX system.
%   CVX_WHERE returns a string containing the base directory of the CVX
%   modeling framework. Within that directory are some useful
%   subdirectories and files:
%       functions/    new functions 
%       examples/     sample cvx models
%       COPYING.txt   copyright information
%   The proper operation of this function assumes that it has not been
%   moved from its default position within the cvx distribution.

s = '-completenames';
s = eval( 'dbstack(s)', 'dbstack' );
s = s(1);
s = eval( 's.file', 's.name' );
if ispc, fs = '\'; else fs = '/'; end
temp = strfind( s, fs );
s( temp(end-1) + 1 : end ) = [];

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
