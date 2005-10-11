function cvx_begin_set

% CVX_BEGIN_SET    Starts a new cvx set specification.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
cvx_setpath( 1 );
cvx_create_problem( true );
assignin( 'caller', 'cvx_problem', cvx_problem );
