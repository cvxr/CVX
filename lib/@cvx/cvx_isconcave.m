function ans = cvx_isconcave( x, full )
error( nargchk( 1, 2, nargin ) );
ans = cvx_vexity( x );
if nargin < 2,
    ans = nnz( ans > 0 ) == 0;
else
    ans = cvx_reshape( ans <= 0, x.size_ );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
