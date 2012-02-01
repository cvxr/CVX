function y = sum_square_abs( x, dim )

%SUM_SQUARE_ABS   Sum of the squares of absolute values.
%   For real arrays, SUM_SQUARE_ABS(X) computes the same result as
%   SUM_SQUARE(X). For complex arrays, SUM_SQUARE(X) first computes the
%   magnitudes of the elements of X, so it compute SUM_SQUARE_ABS(X).
%
%   Similarly, SUM_SQUARE_ABS(X,DIM) implements SUM_SQUARE(ABS(X),DIM).
%
%   Disciplined convex programming information:
%       SUM_SQUARE_ABS(X,...) is convex and nonmonotonic in X. Thus, when
%       used in CVX expressions, X must be affine. DIM must be constant.

error( nargchk( 1, 2, nargin ) );
y = conj( x ) .* x;
if nargin == 2,
    y = sum( y, dim );
else
    y = sum( y );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
