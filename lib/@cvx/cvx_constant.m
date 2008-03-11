function ans = cvx_constant( x )
ans = cvx_reshape( x.basis_( 1, : ), x.size_ );

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
