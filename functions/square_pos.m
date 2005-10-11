function cvx_optval = square_pos( x )
error( nargchk( 1, 1, nargin ) );

%SQUARE    Square of the positive part.
%
%   SQUARE_POS(X) is the square of the postive parts of the elements of X;
%   i.e., SQUARE_POS(X)=MAX(X,0).^2. X must be real.
%
%   Disciplined quadratic programming information:
%       SQUARE_POS(X) is convex and nondecreasing in X. Thus when used in
%       CVX expressions, X must be convex (or affine).

cvx_begin
   variable x2( size( x ) )
   minimize square( x2 )
   x2 >= x;
cvx_end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
