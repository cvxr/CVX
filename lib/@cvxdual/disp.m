function disp( x, prefix, usename )
if nargin < 2, prefix = ''; end
if nargin < 3 | usename,
    disp( [ prefix, 'cvx dual variable ', x.name_, ' (', type( x ), ')' ] );
else,
    disp( [ prefix, 'cvx dual variable (', type( x ), ')' ] );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
