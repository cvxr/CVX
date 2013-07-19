function z = colon( x, y )
if ~isa( x, 'cvxdual' ),
    z = x; x = y; y = z;
end
global cvx___
try
    dvars = cvx___.problems( x.problem_ ).dvars;
    q = subsref( dvars, x.name_ );
catch
    error( 'CVX:Corrupt', 'Internal CVX data corruption. Please CLEAR ALL and rebuild your model.' );
end
if q > 0,
    nm = cvx_subs2str( x.name_ );
    error( 'CVX:DualInUse', 'Dual variable "%s" has already been assigned.', nm(2:end) );
end
dvars = subsasgn( dvars, x.name_, -q );
cvx___.problems( x.problem_ ).dvars = dvars;
try
    z = cvx_setdual( y, x.name_ );
catch
    error( 'CVX:CannotAttachDual', 'Cannot attach a dual variable to an object of type %s.', class( y ) );
end
% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
