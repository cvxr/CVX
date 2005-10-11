function ne( x, y )

%
% Check problem
%

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( cvx_problem, 'cvxprob' ),
    error( 'The ''~='' operator is not defined for structs outside of cvx.' );
end

%
% Always invalid
%

error( sprintf( 'Disciplined convex programming error:\n    Non-equality constraints are not permitted.' ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
