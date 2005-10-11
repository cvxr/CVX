function y = cvx_s_toeplitz( m, n )

if m ~= n,
    error( 'Symmetric structure requires square matrices.' );
end

mn = m * n;
c  = 0 : n - 1;
c  = c( ones( 1, m ), : );
r  = [ 0 : m - 1 ]';
r  = r( :, ones( 1, n ) );
y  = r - c + 1;
temp = y <= 0;
y( temp ) = m + 1 - y( temp );
y  = sparse( 1 : mn, y, 1, mn, m + n - 1 );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.



