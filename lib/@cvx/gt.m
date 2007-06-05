function a = gt( x, y )
error( nargchk( 2, 2, nargin ) );

try
    newcnstr( evalin( 'caller', 'cvx_problem', '[]' ), x, y, '>' );
catch
    error( cvx_lasterr );
end
if nargout > 0,
    a = 'Constraint accepted';
elseif ~isreal( x ) | ~isreal( y ),
    if ~cvx_problem.sdp | ( ( size( x, 1 ) == 1 | size( x, 2 ) == 1 ) & ( size( y, 1 ) == 1 | size( y, 2 ) == 1 ) ),
        error( sprintf( 'Disciplined convex programming error:\n    Both sides of an inequality constraint must be real.' ) );
    end
end    

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
