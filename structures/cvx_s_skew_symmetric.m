function y = cvxs_skew_symmetric( m, n )

if m ~= n,
    error( 'Skew symmetric structure requires square matrices.' );
end

nsq = n * n;
ntr = 0.5 * ( nsq + n );
c  = 0 : n - 1;
c  = c( ones( 1, n ), : );
r  = c';
mn = min( r, c );
mx = max( r, c );
y  = mx + mn .* ( n - 0.5 * ( mn + 1 ) ) + 1;
y  = sparse( 1 : nsq, y( : ), sign( r( : ) - c( : ) ), nsq, ntr );
y  = y( :, any( y, 1 ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
