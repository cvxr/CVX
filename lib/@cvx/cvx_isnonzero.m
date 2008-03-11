function ans = cvx_isnonzero( x, full )
error( nargchk( 1, 2, nargin ) );
ans = any( x.basis_, 1 );
if nargin < 2,
    ans = all( ans );
else
    ans = cvx_reshape( ans, x.size_ );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
