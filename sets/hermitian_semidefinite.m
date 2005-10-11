function cvx_optpnt = hermitian_semidefinite( n )
error( nargchk( 1, 1, nargin ) );
cvx_optpnt = semidefinite( n, true );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
