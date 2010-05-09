function y = cvx_s_skew_symmetric( m, n )
%CVX_S_SKEW_SYMMETRIC Skew-symmetric matrices.

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
y  = sparse( y( : ), 1 : nsq, sign( r( : ) - c( : ) ), ntr, nsq );
y  = y( any( y, 2 ), : );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
