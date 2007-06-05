function y = cvx_s_scaled_identity( m, n )
y = sparse( 1, 1 : m + 1 : m * n, 1, 1, m * n );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
