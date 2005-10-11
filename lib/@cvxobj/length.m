function n = length( x )
error( cvx_verify( x ) );
s = size( x );
if any( s == 0 ),
   n = 0;
else,
   n = max( s );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
