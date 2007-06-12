function v = cvx_vexity( x )

global cvx___
if isempty( x ),
    v = cvx_zeros( x.size_ );
    return
end
p  = cvx___.vexity;
b  = x.basis_;
n  = length( p );
nb = size( b, 1 );
if nb < n,
    p = p( 1 : nb, : );
    b = b( p ~= 0, : );
elseif n < nb,
    p( n, : ) = 0;
    b = b( p ~= 0, : );
else
    b = b( p ~= 0, : );
end
if isempty( b ),
    v = cvx_zeros( x.size_ );
    return
end
p = nonzeros( p ).';
v = p * b;
v( abs( v ) ~= abs( p ) * abs( b ) ) = NaN;
v = sign( v );
v = cvx_reshape( v, x.size_ );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
