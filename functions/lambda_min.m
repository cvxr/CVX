function z = lambda_min( Y )

% LAMBDA_MIN   Minimum eigenvalue of a symmetric matrix.
%     For square matrix X, LAMBDA_MIN(X) is MIN(EIG(X)) if X is Hermitian
%     or symmetric and real; and +Inf otherwise.
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         LAMBDA_MIN is concave and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

error( nargchk( 1, 1, nargin ) );
if ndims( Y ) > 2 | size( Y, 1 ) ~= size( Y, 2 ),
    error( 'Input must be a square matrix.' );
end
err = Y - Y';
Y   = 0.5 * ( Y + Y' );
if norm( err, 'fro' )  > 8 * eps * norm( Y, 'fro' ),
    z = Inf;
else
    z = min( eig( full( Y ) ) );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

