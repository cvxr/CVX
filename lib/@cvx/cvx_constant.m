function ans = cvx_constant( x )
error( cvx_verify( x ) );
ans = cvx_basis( x, 0 );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
