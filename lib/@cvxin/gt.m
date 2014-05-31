function a = gt( x, y )

if isa( y, 'cvxin' ) || ~isa( x, 'cvxin' ) || ~x.active,
	cvx_throw( 'Improper use of the <in> pseudo-operator.' );
end
evalin( 'caller', 'cvx_verify' );
b = cvx_pushcnstr( x.value, y, '==' );
if nargout, a = b; end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
