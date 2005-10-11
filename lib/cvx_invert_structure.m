function x = cvx_invert_structure( x, type )
iscplx = ~isreal( x );
if iscplx,
    x = cvx_c2r( x, 1 );
end
[ m, n ] = size( x );
if nargin == 2 & isequal( type, 'compact' ),
    [ r, c, v ] = find( x );
    temp = [ 1 ; diff( c ) ] ~= 0;
    x = sparse( r( temp ), c( temp ), 1.0 ./ v( temp ), m, n );
else,
    x = x * sparse( 1 : n, 1: n, 1.0 ./ diag( x' * x ), n, n );
end
if iscplx,
    x = cvx_r2c( x, 1 );
end
x = x';

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
