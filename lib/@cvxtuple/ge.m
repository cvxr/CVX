function a = ge( x, y )

evalin( 'caller', 'cvx_verify' );
b = cvx_pushcnstr( x, y, '>=' );
if nargout, a = b; end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
