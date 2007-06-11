function cvx_optval = sum_square( x, dim )
error( nargchk( 1, 2, nargin ) );

%SUM_SQUARE   Sum of squares.
%
%   For vectors, SUM_SQUARE(X) is the sum of the squares of the elements of
%   the vector; i.e., SUM(X.^2).
%
%   For matrices, SUM_SQUARE(X) is a row vector containing the application
%   of SUM_SQUARE to each column. For N-D arrays, the SUM_SQUARE operation
%   is applied to the first non-singleton dimension of X.
%
%   SUM_SQUARE(X,DIM) takes the sum along the dimension DIM of X.
%
%   Disciplined convex programming information:
%       If X is real, then SUM_SQUARE(X,...) is convex and nonmonotonic in
%       X. If X is complex, then SUM_SQUARE(X,...) is neither convex nor
%       concave. Thus, when used in CVX expressions, X must be affine. DIM
%       must be constant.

if ~isreal( x ),
    error( sprintf( 'Disciplined convex programming error:\n   The argument to SUM_SQUARE must be real and affine.' ) );
end
if nargin == 1,
    cvx_optval = quad_over_lin( x, 1 );
else
    cvx_optval = quad_over_lin( x, 1, dim );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
