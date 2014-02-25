function v = cvx_vexity( x )

global cvx___
sz = x.size_;
if any( sz == 0 ),
    v = cvx_zeros( sz );
    return
end
p  = cvx___.vexity;
b  = x.basis_;
n  = length( p );
nb = size( b, 1 );
if nb < n,
    p = p( 1 : nb, 1 );
elseif n < nb,
    p( nb, 1 ) = 0;
end
b = b( p ~= 0, : );
if isempty( b ),
    v = cvx_zeros( sz );
    if x.slow_
        v( isnan( b( 1, : ) ) ) = NaN;
    end
    return
end
s = cvx___.sign;
if nb < n,
    s = s( 1 : nb, 1 );
elseif n < nb,
    s( nb, 1 ) = 0;
end
bs = x.basis_;
ib = ~isreal( bs );
if ib,
    isreal( bs ),
    bi = any( imag( bs ) );
end
s0 = any( bs( s == 0, : ), 1 );
bs = bs( s ~= 0, : );
if cvx___.nan_used,
    bs = sparse( bs );
    b = sparse( b );
end
s = nonzeros(s).';
p = nonzeros(p).';
v = full( p * b );
w = full( s * bs );
w( abs( w ) ~= abs( s ) * abs( bs ) | s0 ) = 0;
tt = abs( v ) ~= abs( p ) * abs( b );
if x.slow_,
    v( tt | isnan( x.basis_( 1, : ) ) ) = NaN;
else
    v( tt ) = NaN;
end
if ib,
    v( bi & v ) = NaN;
    w( bi ) = 0;
end
v = sign( v );
v = v .* ( 1 + ( v == sign( w ) ) );
v = reshape( v, sz );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
