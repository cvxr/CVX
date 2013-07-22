function [ y, symm ] = cvx_s_upper_triangular( cvxm, n, symm )

%CVX_S_UPPER_TRIANGULAR Upper triangular matrices.

[ y, symm ] = cvx_s_banded( m, n, symm, 0, n );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
