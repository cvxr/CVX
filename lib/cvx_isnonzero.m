function y = cvx_isnonzero( x, full ) %#ok
if nargin == 1,
	y = nnz( x ) ~= 0;
else
    y = x ~= 0;
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
