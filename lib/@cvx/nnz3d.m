function s = nnz( x )
error( cvx_verify( x ) );
b = cvx_basis( x );
s = nnz( b ) - nnz( b( :, 1 ) );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
