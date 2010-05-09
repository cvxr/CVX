function y = cvx_s_complex( m, n )
%CVX_S_COMPLEX Complex variables.

mn   = m * n;
rvec = 1 : 2 * mn;
cvec = ceil( 0.5 * rvec );
vvec = ones( 1, 2 * mn );
vvec( 2 : 2 : end ) = 1i;
y = sparse( rvec( : ), cvec( : ), vvec( : ), 2 * mn, mn );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
