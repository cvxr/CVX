function [ y, symm ] = cvx_s_upper_hessenberg( m, n, symm )

%CVX_S_UPPER_HESSENBERG Upper Hessenberg matrices.

[ y, symm ] = cvx_s_banded( m, n, symm, 1, n );

% Copyright 2005-2013 CVX Research, Inc. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
