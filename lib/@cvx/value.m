function v = value( x, data )
error( cvx_verify( x ) );
if nargin == 1,
    data = problem( x );
    data = data.x;
end
s = x.size_;
b = x.basis_;
v = reshape( b * data( 1 : size( b, 2 ), : ), s );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
