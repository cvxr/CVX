function v = vec( x )
error( cvx_verify( x ) );
v = reshape( x, prod( x.size_ ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
