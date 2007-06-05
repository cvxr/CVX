function v = cvx( s, b, d, clean )

switch nargin,
    case 2,
        d = [];
        clean = false;
    case 3,
        clean = true;
    case 1,
        if isa( s, 'cvx' ),
            v = s;
        else
            v = class( struct( 'size_', size( s ), 'basis_', sparse( s(:).' ), 'dual_', '', 'dof_', [] ), 'cvx', cvxobj );
        end
        return
    case 0,
        v = class( struct( 'size_', [ 0, 0 ], 'basis_', sparse( 1, 0 ), 'dual_', '', 'dof_', [] ), 'cvx', cvxobj );
        return
end

if isempty( b ),
    b = sparse( 1, prod( s ) );
elseif issparse( b ) & ~cvx_use_sparse( b ),
    b = full( b );
end
if length( s ) == 1,
    s( 2 ) = 1;
end
if clean,
    b1 = b( 1, : );
    if nnz( b ) == nnz( b1 ),
        v = cvx_reshape( b1, s );
        return
    end
end

v = class( struct( 'size_', s, 'basis_', b, 'dual_', '', 'dof_', d ), 'cvx', cvxobj );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
