function y = sum_square_pos( x, dim )

%SUM_SQUARE_POS   Sum of squares of the positive parts.
%   For vectors, SUM_SQUARE(X) is the sum of the squares of the positive
%   parts of X; i.e., SUM( MAX(X,0)^2 ). X must be real.
%
%   For matrices, SUM_SQUARE_POS(X) is a row vector containing the
%   application of SUM_SQUARE_POS to each column. For N-D arrays, the
%   SUM_SQUARE_POS operation is applied to the first non-singleton
%   dimension of X.
%
%   SUM_SQUARE_POS(X,DIM) takes the sum along the dimension DIM of X.
%
%   Disciplined convex programming information:
%       SUM_SQUARE_POS(X,...) is convex and nondecreasing in X. Thus, when
%       used in CVX expressions, X must be convex (or affine). DIM must
%       always be constant.

error( nargchk( 1, 2, nargin ) );
if nargin == 2,
    y = sum( square_pos( x ), dim );
else
    y = sum( square_pos( x ) );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
