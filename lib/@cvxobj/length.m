function n = length( x )
s = size( x );
if any( s == 0 ),
   n = 0;
else
   n = max( s );
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
