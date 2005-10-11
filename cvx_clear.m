function cvx_clear( arg )

% CVX_CLEAR	Clears all active cvx data in case of an error.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

prob = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( prob, 'cvxprob' ), return; end
p = index( prob );
global cvx___
prob = cvx___.problems( p );
if isempty( prob.variables ),
    varf = {};
else,
    varf = fieldnames( prob.variables );
end
if isempty( prob.duals ),
    dulf = {};
else,
    dulf = fieldnames( prob.duals );
end
evalin( 'caller', sprintf( '%s ', 'clear cvx_problem cvx___ cvx_status cvx_optval', varf{:}, dulf{:} ) );
cvx___.stack( prob.stackpos : end ) = [];
cvx___.problems( p : end ) = [];
if nargin == 0,
    cvx_clearpath( 1 );
end
