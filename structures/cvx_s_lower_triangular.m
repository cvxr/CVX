function [ y, symm ] = cvx_s_lower_triangular( m, n, symm )

%CVX_S_LOWER_TRIANGULAR Lower triangular matrices.

[ y, symm ] = cvx_s_banded( m, n, symm, m, 0 );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
