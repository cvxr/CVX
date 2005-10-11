function y = square_pos( x )
error( nargchk( 1, 1, nargin ) );

%SQUARE    Square of the positive part.
%
%   SQUARE_POS(X) is the square of the postive parts of the elements of X;
%   i.e., SQUARE_POS(X)=MAX(X,0).^2. X must be real.
%
%   Disciplined quadratic programming information:
%       SQUARE_POS(X) is convex and nondecreasing in X. Thus when used in
%       CVX expressions, X must be convex (or affine).

if ~isreal( x ), 
    error( 'Argument must be real.' ); 
end

y = square( max( x, 0 ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
