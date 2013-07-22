function [ y, symm ] = cvx_s_toeplitz( m, n, symm )
%CVX_S_TOEPLITZ Toeplitz matrices.
mn = m * n;
c  = 0 : n - 1;
c  = c( ones( 1, m ), : );
r  = ( 0 : m - 1 )';
r  = r( :, ones( 1, n ) );
y  = r - c + 1;
temp = y <= 0;
if symm,
    y( temp ) = 2 - y(temp);
    maxy = max( m, n );
    symm = false;
else
    y( temp ) = m + 1 - y( temp );
    maxy = m + n - 1;
end
y = sparse( y, 1 : mn, 1, maxy, mn );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
