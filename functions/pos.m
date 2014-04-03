function y = pos( x )

% POS    Positive part.
%    POS(X) = MAX(X,0). X must be real.
%
%    Disciplined convex programming information:
%        POS(X) is convex and nondecreasing in X. Thus when used in CVX
%        expressions, X must be convex (or affine).

if ~isreal( x ), error( 'Argument must be real.' ); end
y = max( x, 0 );

% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
