function z = isequal( x, y )
error( cvx_verify( x ) );
z = x.id_ == y.id_;
if ~z,
    error( cvx_verify( y ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
