function y = inv_pos( x )

%INV_POS   Reciprocal of a positive quantity.
%    INV_POS(X) returns 1./X if X is positive, and +Inf otherwise.
%    X must be real.
%
%    For matrices and N-D arrays, the function is applied to each element.
%
%     Disciplined convex programming information:
%         INV_POS is convex and nonincreasing; therefore, when used in CVX
%         specifications, its argument must be concave (or affine).

y = recip( pdom( x ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
