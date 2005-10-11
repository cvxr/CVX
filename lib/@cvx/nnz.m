function s = nnz( x )
error( cvx_verify( x ) );
s = nnz( any( cvx_basis( x ), 2 ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
