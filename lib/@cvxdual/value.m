function v = value( x )
v = cvxaff( x );
if isa( v, 'cvx' ),
    global cvx___
    v = value( v, cvx___.y );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
