function cvx_optval = pos( x )
error( nargchk( 1, 1, nargin ) );

% POS    Positive part.
%    
%    POS(X) = MAX(X,0).
%
%    Disciplined convex programming information:
%        POS(X) is convex and nondecreasing in X. Thus when used in CVX
%        expressions, X must be convex (or affine).

cvx_optval = max( x, 0 );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
