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
s0 = any( b( s == 0, : ), 1 ) | any( imag( b ), 1 );
b  = b( s ~= 0, : );
s  = nonzeros(s).';
w  = full( s * b );
w( abs( w ) ~= abs( s ) * abs( b ) | s0 ) = 0;
if nargin == 2, 
    w( ~any( x.basis_, 1 ) ) = dz;
end
w = reshape( sign( w ), sz );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
