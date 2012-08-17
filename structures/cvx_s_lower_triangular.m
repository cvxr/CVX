function y = cvx_s_lower_triangular( m, n )
%CVX_S_LOWER_TRIANGULAR Lower triangular matrices.
y = cvx_s_banded( m, n, m, 0 );

% Copyright 2012 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
