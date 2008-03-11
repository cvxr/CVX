function z = lambda_min( x )

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
if ndims( x ) > 2 | size( x, 1 ) ~= size( x, 2 ),

    error( 'Input must be a square matrix.' );

elseif any( any( x ~= x' ) ),

    z = Inf;

else

    z = min( eig( full( x ) ) );

end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

