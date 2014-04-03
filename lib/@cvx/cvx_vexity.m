function [ v, w ] = cvx_vexity( x )

% Vexity:
%  -1: concave
%   0: affine
%  +1: convex
% NaN: invalid

global cvx___
sz = x.size_;
if any( sz == 0 ),
    v = cvx_zeros( sz );
    w = v;
    return
end
p  = cvx___.vexity;
b  = x.basis_;
ib = ~isreal( b );
if ib, bi = any( imag( b ), 1 ); end
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
else
    if cvx___.nan_used, b = sparse( b ); end
    p = nonzeros(p).';
    v = full( p * b );
    tt = abs( v ) ~= abs( p ) * abs( b );
    if ib
        tt = tt | ( ib & v );
    end
    if x.slow_
        tt = tt | isnan( x.basis_( 1, : ) );
    end
    v( tt ) = NaN;
    v = reshape( sign( v ), sz );
end
if nargout < 2,
    return
end

%
% New addition: compute the sign 
%  -1: nonpositive
%   0: unknown sign
%  +1: nonnegative
%

s = cvx___.sign;
if nb < n,
    s = s( 1 : nb, 1 );
elseif n < nb,
    s( nb, 1 ) = 0;
end
b = x.basis_;
tt = any( b( s == 0, : ), 1 );
b = b( s ~= 0, : );
if ~isempty( b ),
    s = nonzeros(s).';
    w = full( s * b );
    tt = tt | abs( w ) ~= abs( s ) * abs( b );
end
if ib, tt = tt | bi; end
w( tt ) = 0;
w = reshape( w, sz );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
