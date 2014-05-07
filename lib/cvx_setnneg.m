function cvx_setnneg( tx )

global cvx___
if isa( tx, 'cvx' ),
    tx = find( any( cvx_basis( tx ), 2 ) );
end
persistent mpos;
if isempty( mpos ),
	mpos = int8([2,2,3,22,6,7,7,9,10,10,12,13,13,22,15,16,17,18,19,20,21,22]);
end
cvx___.classes(tx) = mpos(cvx___.classes(tx));

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
