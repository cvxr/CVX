function z = isequal( x, y )

% This only works if we treat CVX objects as immutable.
z = x.id_ == y.id_;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

