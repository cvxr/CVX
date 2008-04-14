function cvx_optval = lambda_max( x )

%LAMBDA_MIN   Internal cvx version.

error( nargchk( 1, 1, nargin ) );
if ndims( x ) > 2 | size( x, 1 ) ~= size( x, 2 ),

    error( 'Input must be a square matrix.' );

elseif cvx_isconstant( x ),

    cvx_optval = cvx( lambda_min( cvx_constant( x ) ) );

elseif cvx_isaffine( x ),

    n = size( x, 1 );
    cvx_begin
        hypograph variable z
        x - z * eye( n ) == semidefinite( n, ~isreal( x ) );
    cvx_end

else

    error( 'Discipliend convex programming error:\n    LAMBDA_MIN is concave and nonmonotonic, so its input must be affine.' );

end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
