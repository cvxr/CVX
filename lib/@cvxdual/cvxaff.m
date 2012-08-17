function y = cvxaff( x )
global cvx___
y = cvx___.problems( index( x.problem_ ) );
y = builtin( 'subsref', y.duals, x.name_ );

% Copyright 2012 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
