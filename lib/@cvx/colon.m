function z = colon( x, y )

if ~isa( x, 'cvxdual' ),
    z = cvx( x );
    d = y;
elseif ~isa( y, 'cvxdual' ),
    z = cvx( y );
    d = x;
else
    error( 'Usage: <dual var> : <expression> or <expression> : <dual var>' );
end

if inuse( d ),
    nm = cvx_subs2str( name( d ) );
    error( [ 'Dual variable "', nm(2:end), '" has already been used.' ] );
else
    z = setdual( z, name( d ) );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
