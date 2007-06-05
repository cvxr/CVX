function ans = cvx_isconvex( x, full )
error( nargchk( 1, 2, nargin ) );
if nargin < 2,
    ans = nnz( ans < 0 ) == 0;
else
    ans = cvx_reshape( ans >= 0, x.size_ );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
