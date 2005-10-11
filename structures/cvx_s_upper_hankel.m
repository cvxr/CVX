function y = cvx_s_hankel( m, n )
c  = 0 : n - 1;
c  = c( ones( 1, m ), : );
r  = [ 0 : m - 1 ]';
r  = r( :, ones( 1, n ) );
v  = abs( r + c ) + 1;
temp = v <= min( m, n );
y = sparse( r( temp ) + m * c( temp ) + 1, v( temp ), 1, m * n, min( m, n ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
