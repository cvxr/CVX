function x = vec( x )

% VEC   CVX implementation of vec

s = x.size_;
n = prod(s);
if s(1) ~= n,
	x = cvx( [ n, 1 ], x.basis_ );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
