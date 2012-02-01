function v = cvxcnst( p, rhs )
v = class( struct( 'problem', p, 'rhs', rhs ), 'cvxcnst' );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
