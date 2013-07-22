function [ y, symm ] = cvx_s_tridiagonal( m, n, symm )

%CVX_S_TRIDIAGONAL Tridiagional matrices.

[ y, symm ] = cvx_s_banded( m, n, symm, 1, 1 );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
