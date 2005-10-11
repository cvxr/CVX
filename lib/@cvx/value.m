function v = value( x, data )
error( cvx_verify( x ) );
if nargin == 1,
    data = problem( x );
    data = data.x;
end
s = x.size_;
b = x.basis_;
v = reshape( b * data( 1 : size( b, 2 ), : ), s );
if ~isempty( x.geom_ ),
    g = x.geom_;
    l = prod( s );
    g = reshape( g * data( 1 : size( g, 2 ), : ), l, size( g, 1 ) / l );
    v = v + reshape( exp( sum( g, 2 ) ), s );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
