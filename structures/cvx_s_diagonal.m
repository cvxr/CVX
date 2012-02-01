function y = cvx_s_diagonal( m, n )
%CVX_S_DIAGONAL Diagonal matrices.
y = cvx_s_banded( m, n, 0, 0 );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
