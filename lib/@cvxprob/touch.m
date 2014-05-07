function x = touch( prob, x, iseq ) %#ok
if ~isa( x, 'cvx' ), return; end
global cvx___
p = prob.index_;
pstr = cvx___.problems( p );
v = pstr.t_variable;
nv = numel( v );
if nv <= 1, return; end
b  = cvx_basis( x );
y  = any( b, 2 );
ny = numel( y );
if ny < nv,
    v( find( y ) ) = 1; %#ok
else
    v( y( 1 : nv ) ) = 1;
end
cvx___.problems( p ).t_variable = v;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
