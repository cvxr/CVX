function y = pow_abs( x, p )

%POW_ABS   Internal cvx version.

error( nargchk( 2, 2, nargin ) );
if ~isa( p, 'double' ) | ~isreal( p ) | nnz( isinf( p ) | isnan( p ) ) | nnz( p < 1 ),
    error( 'Second argument must be greater than or equal to 1.' );
end

y = pow_pos( abs( x ), p );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

