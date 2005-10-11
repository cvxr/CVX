function ans = cvx_isconstant( x, full )
error( nargchk( 1, 2, nargin ) );
b = cvx_basis( x );
if nargin < 2,
    ans = nnz( b ) == nnz( b( :, 1 ) );
else,
    ans = reshape( ~any( b( :, 2 : end ), 2 ), size( x ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
