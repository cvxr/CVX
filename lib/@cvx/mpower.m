function z = mpower( x, y )

if length( x ) == 1 & length( y ) == 1,
    z = power( x, y );
else
    error( sprintf( 'Disciplined convex programming error:\n    Matrix powers not permitted.' ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
