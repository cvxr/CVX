function z = colon( x, y )

if ~isa( x, 'cvxdual' ),
    error( cvx_verify( x ) );
    z = cvx( x );
    d = y;
elseif ~isa( y, 'cvxdual' ),
    error( cvx_verify( y ) );
    z = cvx( y );
    d = x;
else,
    error( 'Usage: <dual var> : <expression> or <expression> : <dual var>' );
end

if inuse( d ),
    error( [ 'Dual variable "', name( d ), '" has already been used.' ] );
else,
    z = setdual( z, name( d ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
