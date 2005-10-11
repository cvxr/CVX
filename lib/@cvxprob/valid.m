function s = valid( x )

global cvx___
ndx = index( x );
if ndx > length( cvx___.problems ),
    s = false;
else,
    s = id( x ) == id( cvx___.problems( ndx ).self );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
