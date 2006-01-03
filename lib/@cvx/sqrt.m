function cvx_optval = sqrt( x )
error( nargchk( 1, 1, nargin ) );

%SQRT    Square root.
%
%   SQRT(X) is the square root of the elements of X. If X is numeric, then
%   complex results are returned for those elements of X which are complex
%   or real and negative. If X is a CVX expression, then X must be real,
%   the function effectively constraints X to be nonnegative.
%
%   Disciplined quadratic programming information:
%       SQRT(X) is concave and nondecreasing in X. Thus when used in CVX
%       expressions, X must be concave (or affine).

sx = size( x );
cvx_begin
    variable y( sx )
    maximize y
    square( y ) <= x;
cvx_end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
