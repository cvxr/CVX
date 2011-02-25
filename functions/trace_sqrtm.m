function z = trace_sqrtm( Y )

% TRACE_SQRTM   Trace of the square root of a PSD matrix.
%     For square matrix X, TRACE_SQRTM(X) is TRACE(SQRTM(X)) if X is Hermitian
%     or symmetric and positive definite; and +Inf otherwise.
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         TRACE_SQRTM is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

error( nargchk( 1, 1, nargin ) );
if ndims( Y ) > 2 || size( Y, 1 ) ~= size( Y, 2 ),
    error( 'Input must be a square matrix.' );
end
err = Y - Y';
Y   = 0.5 * ( Y + Y' );
if norm( err, 'fro' )  > 8 * eps * norm( Y, 'fro' ),
    z = Inf;
else
    z = eig( full( Y ) );
    if any( z <= 0 ),
        z = Inf;
    else
        z = sum(sqrt(z));
    end
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

