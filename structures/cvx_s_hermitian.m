function y = cvx_s_hermitian( m, n )

if m ~= n,
    error( 'Hermitian structure requires square matrices.' );
end

nsq = n * n;
ntr = 0.5 * ( nsq + n );
c  = 0 : n - 1;
c  = c( ones( 1, 2 * n ), : );
c  = c( : );
r  = 0 : n - 1;
r  = r( [ 1, 1 ], : );
r  = r( : );
r  = r( :, ones( 1, n ) );
r  = r( : );
v  = [ 1 ; j ];
v  = v( :, ones( 1, nsq ) );
v  = v( : );
temp = r < c;
v( temp ) = conj( v( temp ) );
temp = r == c;
v( temp ) = real( v( temp ) );
mn = min( r, c );
mx = max( r, c );
y  = sparse( r + n * c + 1, 2 * ( mx + mn .* ( n - 0.5 * ( mn + 1 ) ) + 1 ) - ( v == 1 ), v, nsq, length( v ) );
y  = y( :, any( y, 1 ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
