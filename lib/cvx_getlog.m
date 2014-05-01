function by = cvx_getlog( bx )

% WARNING: This assumes that there is exactly one non-zero element per
% column, that all non-zeros are positive, and that any element appearing
% in rows greater than 1 has a non-zero entry in cvx___.logarithm.

global cvx___
nb = size( bx, 2 );
[ rx, cx, vx ] = find( bx );
tt = rx > 1;
ry = rx(tt);
cy = cx(tt);
logs = cvx___.logarithm( ry, 1 );
nl = max([1,max(logs)]);
by = sparse( logs, cy, 1, nl, nb ) + sparse( 1, cx, log(vx), nl, nb );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
