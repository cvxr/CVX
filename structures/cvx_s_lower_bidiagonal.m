function y = cvx_s_lower_bidiagonal( m, n )
%CVX_S_LOWER_BIDIAGONAL Lower bidiagonal matrices.
y = cvx_s_banded( m, n, 1, 0 );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
