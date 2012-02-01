function y = dof( x )
y = size( cvx_basis( cvxaff( x ) ), 2 );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
