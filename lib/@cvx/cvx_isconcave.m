function ans = cvx_isconcave( x, full )
error( nargchk( 1, 2, nargin ) );
error( cvx_verify( x ) );
ans = cvx_vexity( x ) <= 0;
if nargin < 2,
    ans = all( ans( : ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
