function v = cvx( p, s, b, d )
switch nargin,
    case 0,
        v = cvx( cvxprob( 'current' ), [] );
    case 1,
        switch class( p ),
            case 'cvx',
                error( cvx_verify( p ) );
                v = p;
            case 'cvxprob',
                v = cvx( p, [] );
            case 'double',
                v = cvx( cvxprob( 'current' ), p );
            otherwise,
                error( sprintf( 'Cannot create a cvx object from class %s', class( p ) ) );
        end
    case 2,
        if ~isa( p, 'cvxprob' ) | length( p ) ~= 1,
            error( 'First argument must be a scalar cvx problem object.' );
        elseif ~isa( s, 'double' ) & ~isa( s, 'sparse' ),
            error( 'Second argument must be numeric.' );
        elseif any( isinf( s( : ) ) ) | any( isnan( s( : ) ) ),
            error( 'Infs and NaNs are not permitted.' );
        end
        v = class( struct( 'size_', size( s ), 'basis_', s( : ), 'dual_', '', 'dof_', [] ), 'cvx', cvxobj( p ) );
    case { 3, 4 },
        if ~isa( p, 'cvxprob' ),
            error( 'First argument must be a scalar cvx problem object.' );
        end
        [ temp, s ] = cvx_check_dimlist( s, true );
        if ~temp,
            error( 'Second argument must be a dimension list.' );
        elseif ~isa( b, 'double' ) | ndims( b ) > 2 | size( b, 1 ) ~= prod( s ),
            error( 'Third argument must be a basis matrix.' );
        elseif nargin < 4 | isempty( d ),
            d = [];
        elseif ~isa( d, 'double' ) | length( d ) ~= 1 | ~isreal( d ) | d < 0,
            error( 'Fourth argument must be an integer.' );
        end
        if isempty( b ),
            b = sparse( [], [], [], prod( s ), 1 );
        elseif ~issparse( b ) & nnz( b ) * 4 < prod( size( b ) ),
            b = sparse( b );
        end
        v = class( struct( 'size_', s, 'basis_', b, 'dual_', '', 'dof_', d ), 'cvx', cvxobj( p ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
