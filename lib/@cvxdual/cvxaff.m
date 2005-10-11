function y = cvxaff( x )
y = struct( problem( x ) );
y = subsref( y.duals, struct( 'type', '.', 'subs', x.name_ ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
