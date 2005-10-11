function cvx_optval = sum_square( x, dim )
error( nargchk( 1, 2, nargin ) );

%SUM_SQUARE   Sum of squares.
%
%   For real vectors, SUM_SQUARE(X) is the sum of the squares of the
%   elements of the vector; i.e., SUM( X .^ 2 ).
%
%   For complex vectors, SUM_SQUARE(X) is the sum of the squares of the
%   magnitudes of the vector; i.e., SUM( ABS( X ) .^ 2 ).
%
%   For matrices, SUM_SQUARE(X) is a row vector containing the application
%   of SUM_SQUARE to each column. For N-D arrays, the SUM_SQUARE operation
%   is applied to the first non-singleton dimension of X.
%
%   SUM_SQUARE(X,DIM) takes the sum along the dimension DIM of X.
%
%   Disciplined convex programming information:
%       SUM_SQUARE(X,...) is convex and nonmonotonic in X. Thus, when used
%       in CVX expressions, X must be affine. DIM must always be constant.

sx = size( x );
if nargin < 2,
    dim = cvx_default_dimension( sx );
end

cvx_begin
    variable x2( sx )
    minimize sum_square( x2, dim )
    x2 >= x;
cvx_end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
