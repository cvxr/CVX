function z = matrix_frac( x,Y )
error( nargchk( 2, 2, nargin ) );

% MATRIX_FRAC
%     For square matrix Y and vector x of same size, MATRIX_FRAC(x,Y)
%     is x'*(inv(Y)*x) if Y is Hermitian
%     or symmetric and real positive definite; +Inf otherwise.
%
%     An error results if Y is not a square matrix, or the size of
%     the vector x does not match the size of matrix Y
%
%     Disciplined convex programming information:
%         MATRIX_FRAC is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

if ndims( Y ) > 2 | size( Y, 1 ) ~= size( Y, 2 ),

    error( 'Second argument must be a square matrix.' );

elseif any( any( Y ~= Y' ) ) | min(eig(Y))<=0,

    z = Inf;

else

    z = x'*(Y\x);

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
