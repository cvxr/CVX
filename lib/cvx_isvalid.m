function y = cvx_isvalid( x, full )
y = isnan( x ) | isinf( x );
if nargin < 2 || ~full,
	y = ~any( y(:) );
else
	y = ~y;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
