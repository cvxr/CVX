function y = pow_abs( x, p )

% POW_ABS   Absolute value raised to a power.
%    POW_POS(X,P) = ABS(X).^P. The dimensions of X and P must be identical,
%    or one must be a scalar. P must be greater than or equal to 1.
%    
%    POW_ABS actually converts P to its nearest rational representation, as
%    given by the RAT() function. Please see the user's guide in the
%    section "Powers, p-norms, and polynomials" for more information about
%    how this function is implemented.
%
%    Disciplined convex programming information:
%        POW_ABS is convex and nonmonotonic in X; therefore, when used in
%        CVX models, X must be affine.

%
% Check second argument to make sure it is >= 1
%

error( nargchk( 2, 2, nargin ) );
if ~isa( p, 'double' ) | ~isreal( p ) | nnz( isinf( p ) | isnan( p ) ) | nnz( p < 1 ),
    error( 'Second argument must be greater than or equal to 1.' );
end

y = pow_pos( abs( x ), p );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

