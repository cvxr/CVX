function cvx_optval = lambda_max( x )
error( nargchk( 1, 1, nargin ) );

% LAMBDA_MIN   Minimum eigenvalue of a symmetric matrix.
%     For square matrix X, LAMBDA_MIN(X) is MIN(EIG(X)) if X is Hermitian
%     or symmetric and real; and +Inf otherwise.
%
%     An error results if X is not a square matrix.
%
%     Disciplined convex programming information:
%         LAMBDA_MIN is concave and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

if ndims( x ) > 2 | size( x, 1 ) ~= size( x, 2 ),
    
    error( 'Input must be a square matrix.' );
    
elseif cvx_constant( x ),
    
    cvx_optval = lambda_min( cvx_constant( x ) );
    
elseif cvx_isaffine( x ),
    
    n = size( x, 1 );
    cvx_begin
        variable z
        maximize z
        subject to
            x - z * eye( n ) == semidefinite( n, ~isreal( x ) );
    cvx_end

else,
    
    error( 'Discipliend convex programming error:\n    LAMBDA_MIN is concave and nonmonotonic, so its input must be affine.' );
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
