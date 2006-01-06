function ans = cvx_isnonzero( x, full )
error( nargchk( 1, 2, nargin ) );
error( cvx_verify( x ) );
ans = any( cvx_basis( x ), 2 );
if nargin < 2,
    ans = all( ans );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
