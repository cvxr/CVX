function y = cvxs_symmetric( m, n )

if m ~= n,
    error( 'Symmetric structure requires square matrices.' );
end

nsq = n * n;
ntr = 0.5 * ( nsq + n );
c  = 0 : n - 1;
c  = c( ones( 1, n ), : );
r  = c';
mn = min( r, c );
mx = max( r, c );
y  = mx + mn .* ( n - 0.5 * ( mn + 1 ) ) + 1;
y  = sparse( 1 : nsq, y( : ), 1, nsq, ntr );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

