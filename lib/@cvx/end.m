function y = end( x, k, n )
error( cvx_verify( x ) );
sz = size( x );
nz = length( sz );
if k > nz,
    y = 1;
elseif k < n | nz <= n,
    y = sz( k );
else,
    y = prod( sz( k : end ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
