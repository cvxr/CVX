function cvx_begin

% CVX_BEGIN    Starts a new cvx specification.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
cvx_setpath( 1 );
cvx_create_problem( false );
assignin( 'caller', 'cvx_problem', cvx_problem );
