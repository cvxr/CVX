function z = matrix_frac( x,Y )

%MATRIX_FRAC   Matrix fractional function.
%     MATRIX_FRAC(x,Y), where Y is a square matrix and x is a vector of the
%     same size, computes x'*(inv(Y)*x) if Y is Hermitian positive definite, and
%     +Inf otherwise.
%
%     An error results if Y is not a square matrix, or the size of
%     the vector x does not match the size of matrix Y.
%
%     Disciplined convex programming information:
%         MATRIX_FRAC is convex and nonmonotonic (at least with respect to
%         elementwise comparison), so its argument must be affine.

error( nargchk( 2, 2, nargin ) );
if ndims( Y ) > 2 | size( Y, 1 ) ~= size( Y, 2 ),

    error( 'Second argument must be a square matrix.' );

elseif any( any( Y ~= Y' ) ) | min(eig(Y))<=0,

    z = Inf;

else

    z = x'*(Y\x);

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
