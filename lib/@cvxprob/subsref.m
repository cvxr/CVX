function x = subsref( x, S, cheat )

try
    global cvx___
    x = subsref( cvx___.problems( index( x ) ), S );
catch
    error( lasterr );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
