function y = square_abs( x )

%SQUARE_ABS   Square of absolute value.
%   For real X, SQUARE_ABS(X) produces the same result as SQUARE(X); that
%   is, it squares the elements of X. For complex X, SQUARE_ABS(X) returns
%   a real array whose elements are the squares of the magnitudes of the
%   elements of X; that is, SQUARE_ABS(X) = CONJ(X).*X.
%
%   Disciplined convex programming information:
%       SQUARE_ABS(X) is convex and nonmonotonic in X. Thus when used in 
%       CVX expressions, X must be affine.

error( nargchk( 1, 1, nargin ) );
y = conj( x ) .* x;

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
