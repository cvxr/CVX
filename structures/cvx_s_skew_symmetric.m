function [ y, symm ] = cvx_s_skew_symmetric( m, n, symm )
%CVX_S_SKEW_SYMMETRIC Skew-symmetric matrices.

if m ~= n,
    error( 'Skew symmetric structure requires square matrices.' );
elseif symm,
	error( 'Cannot specify both symmetry and skew symmetry.' );
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

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
