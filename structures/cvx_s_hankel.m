function y = cvx_s_hankel( m, n )
%CVX_S_HANKEL Hankel matrices.
c  = 0 : n - 1;
c  = c( ones( 1, m ), : );
r  = ( 0 : m - 1 )';
r  = r( :, ones( 1, n ) );
v  = abs( r + c ) + 1;
y = sparse( v, r + m * c + 1, 1, m + n + 1, m * n );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
