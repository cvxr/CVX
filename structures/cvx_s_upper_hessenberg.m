function y = cvx_s_upper_hessenberg( m, n )
%CVX_S_UPPER_HESSENBERG Upper Hessenberg matrices.
y = cvx_s_banded( m, n, 1, n );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
