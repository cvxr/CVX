function x = cvx_invert_structure( x, compact )
iscplx = ~isreal( x );
if iscplx,
    x = cvx_c2r( x, 2 );
end
[ m, n ] = size( x );
if nargin == 2,
    [ c, r, v ] = find( x' );
    temp = [ 1 ; diff( r ) ] ~= 0;
    x = sparse( r( temp ), c( temp ), 1.0 ./ v( temp ), m, n );
else
    x = sparse( 1 : m, 1 : m, 1.0 ./ diag( x * x' ), m, m ) * x;
end
if iscplx,
    x = cvx_r2c( x, 2 );
end
x = x';

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
