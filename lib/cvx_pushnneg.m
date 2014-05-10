function ndim = cvx_pushnneg( nx )

if nx == 0,
    ndim = zeros(0,1);
    return
end
global cvx___
persistent aff_code
if isempty( aff_code ),
    aff_code = int8(find(cvx_remap('p_affine_')));
end
nO = length( cvx___.classes );
ndim = nO + 1 : nO + nx;
cvx___.classes ( ndim, 1 ) = aff_code;
cones = cvx___.cones;
if ~isempty( cones ) && isequal( cones(1).type, 'nonnegative' ),
    cones(1).indices = [ cones(1).indices, ndim ];
else
    ncone = struct( 'type', 'nonnegative', 'indices', ndim, 'slacks', 1 );
    cones = [ ncone, cones ];
end
cvx___.cones = cones;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

