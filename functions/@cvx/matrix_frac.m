function cvx_optval = matrix_frac( x,Y )

%MATRIX_FRAC   Internal cvx version.

error( nargchk( 2, 2, nargin ) );
if ndims( Y ) > 2 || size( Y, 1 ) ~= size( Y, 2 ),

    error( 'Second argument must be a square matrix.' );

elseif ndims( x ) > 2 || size( x, 2 ) > 1,

    error( 'First argument must be a column vector.' );

elseif size( x, 1 ) ~= size( Y, 1 ),

    error( 'Size of first argument (vector) must match size of second argument (matrix).' );

elseif cvx_isconstant( x ) && cvx_isconstant( Y ),

    cvx_optval = cvx( matrix_frac( cvx_constant( x ), cvx_constant(Y) ) );

elseif cvx_isaffine( x ) && cvx_isaffine( Y ),

    n = size( x, 1 );
    cvx_begin
        epigraph variable z
        [Y x; x' z] == semidefinite( n+1 );
    cvx_end

else

    error( 'Disciplined convex programming error:\n    MATRIX_FRAC is convex and nonmonotonic, so its input must be affine.' );

end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
