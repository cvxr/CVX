function y = cvx_s_complex( m, n )

mn = m * n;
rvec = 1 : mn;
rvec = rvec( [ 1, 1 ] , : );
rvec = rvec( : )';
cvec = 1 : 2 * mn;
vvec = ones( 1, 2 * mn );
vvec( 2 : 2 : end ) = j;
y = sparse( rvec( : ), cvec( : ), vvec( : ), mn, 2 * mn );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

