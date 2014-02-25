function w = cvx_sign( x, dz )

global cvx___
sz = x.size_;
if any( sz == 0 ),
    w = cvx_zeros( sz );
    return
end
s  = cvx___.sign;
b  = x.basis_;
n  = length( s );
nb = size( b, 1 );
if nb < n,
    s = s( 1 : nb, 1 );
elseif n < nb,
    s( nb, 1 ) = 0;
end
sn = ~any( b, 1 );
s0 = any( b( ~s, : ) ) | any( imag( b ) );
b  = b( s ~= 0, : );
s  = nonzeros(s).';
if cvx___.nan_used,
    b = sparse( b );
end
w = full( s * b );
w( abs( w ) ~= abs( s ) * abs( b ) | s0 ) = 0;
if nargin == 2, w( sn ) = dz; end
w = reshape( sign( w ), sz );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
