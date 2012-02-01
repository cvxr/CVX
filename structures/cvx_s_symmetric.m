function y = cvx_s_symmetric( m, n )
%CVX_S_SYMMETRIC Symmetric matrices (lower triangle storage).

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
y  = sparse( y( : ), 1 : nsq, 1, ntr, nsq );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

