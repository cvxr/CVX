function z = colon( x, y )
if ~isa( x, 'cvxdual' ),
    z = x; x = y; y = z;
end
global cvx___
try
    duals = cvx___.problems( x.problem_ ).duals;
    q = builtin( 'subsref', duals, x.name_ );
catch
    cvx_throw( 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
end
if isa( q, 'cvx' ),
    nm = cvx_subs2str( x.name_ );
    cvx_throw( 'Dual variable "%s" has already been assigned.', nm(2:end) );
end
try
    z = cvx_setdual( y, x.name_ );
catch
    cvx_throw( 'Cannot attach a dual variable to an object of type %s.', class( y ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
