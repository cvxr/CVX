function v = cvx_classify( x )

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - nonpositive concave
% 6  - concave
% 7  - real affine
% 8  - convex
% 9  - nonnegative convex
% 10 - complex affine
% 11 - log concave
% 12 - log affine
% 13 - log convex monomial
% 14 - log convex posynomial
% 15 - invalid

v = full( sign( real( x ) ) ) + 2;
if ~isreal( x ),
	v( imag( x ) ~= 0 ) = 4;
end
v( ~isfinite( x ) ) = 15;
v = reshape( v, 1, numel( x ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
