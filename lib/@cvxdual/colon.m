function z = colon( x, y )

if ~isa( x, 'cvxdual' ),
    error( cvx_verify( x ) );
    x = cvx_collapse( x );
    switch class( x ),
        case { 'cell', 'struct' },
            error( 'Cannot assign a dual variable to a composite constraint.' );
        otherwise,
            x = cvx( x );
    end
    z = x;
    d = y;
elseif ~isa( y, 'cvxdual' ),
    error( cvx_verify( y ) );
    y = cvx_collapse( y );
    switch class( y ),
        case { 'cell', 'struct' },
            error( 'Cannot assign a dual variable to a composite constraint.' );
        otherwise,
            y = cvx( y );
    end
    z = y;
    d = x;
else,
    error( 'Usage: <dual var> : <expression> or <expression> : <dual var>' );
end

if inuse( d ),
    error( [ 'Dual variable "', d.name_, '" has already been used.' ] );
else,
    z = setdual( z, d.name_ );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
