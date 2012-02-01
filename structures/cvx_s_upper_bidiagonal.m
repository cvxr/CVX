function y = cvx_s_upper_bidiagonal( m, n )
%CVX_S_UPPER_BIDIAGONAL Upper bidiagonal matrices.
y = cvx_s_banded( m, n, 0, 1 );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
