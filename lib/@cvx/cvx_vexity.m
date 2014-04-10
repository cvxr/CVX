function v = cvx_vexity( x )

% Vexity:
%  -1: concave
%   0: affine
%  +1: convex
% NaN: invalid

global cvx___
if 1,
    v = cvx_vexity_mex( x.basis_, cvx___.vexity );
    v = reshape( v, x.size_ );
    return
else
    sz = x.size_;
    if any( sz == 0 ),
        v = cvx_zeros( sz );
        w = v;
        return
    end
    p  = cvx___.vexity;
    b  = x.basis_;
    ib = ~isreal( b );
    if ib, 
        bi = ~any( imag( b ), 1 ); 
    end
    n  = length( p );
    nb = size( b, 1 );
    if nb < n,
        p = p( 1 : nb, 1 );
    elseif n < nb,
        p( nb, 1 ) = 0;
    end
    b = b( p ~= 0, : );
    if isempty( b ),
        v = zeros( sz );
        if x.slow_
            v( isnan( x.basis_( 1, : ) ) ) = NaN;
        end
    else
        if cvx___.nan_used, b = sparse( b ); end
        p = nonzeros(p).';
        v = full( p * b );
        tt = abs( v ) ~= abs( p ) * abs( b );
        if ib
            tt = tt & ( bi | ( v == 0 ) );
        end
        if x.slow_
            tt = tt | isnan( x.basis_( 1, : ) );
        end
        v( tt ) = NaN;
        v = reshape( sign( v ), sz );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
