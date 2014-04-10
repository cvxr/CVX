function cvx_optval = matrix_frac( x, Y )

%MATRIX_FRAC   Internal cvx version.

if ndims( Y ) > 2 || size( Y, 1 ) ~= size( Y, 2 ), %#ok

    error( 'Second argument must be a square matrix.' );

elseif ndims( x ) > 2 || size( x, 2 ) > 1, %#ok

    error( 'First argument must be a column vector or matrix.' );

elseif size( x, 1 ) ~= size( Y, 1 ),

    error( 'The number of rows in X and Y must match.' );

elseif cvx_isconstant( x ) && cvx_isconstant( Y ),

    cvx_optval = cvx( matrix_frac( cvx_constant( x ), cvx_constant(Y) ) );

elseif cvx_isaffine( x ) && cvx_isaffine( Y ),

    [n,m] = size(x);
    z = [];
    cvx_begin
        variable z(n)
        minimize(sum(n))
        [Y x; x' diag(z)] == semidefinite( n+m ); %#ok
        cvx_setnneg( z );
    cvx_end

else

    error( 'Disciplined convex programming error:\n    MATRIX_FRAC is convex and nonmonotonic, so its input must be affine.' );

end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
