function v = cvx_classify( x )

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - negative concave
% 6  - concave
% 7  - positive concave
% 8  - negative affine
% 9  - real affine
% 10 - positive affine
% 11 - negative convex
% 12 - convex
% 13 - positive convex
% 14 - complex affine
% 15 - log concave
% 16 - log affine
% 17 - log convex monomial
% 18 - log convex posynomial
% 19 - invalid

global cvx___
[ v, w ] = cvx_vexity( x );
v = reshape( full( 3 * v + w ), 1, prod( x.size_ ) );
if isempty( x ), return; end
b = x.basis_ ~= 0;
q = sum( b, 1 );
s = b( 1, : );

% Constants
tt = q == s;
if any( tt ),
    if ~isreal( x.basis_ ),
        ti = any( imag( x.basis_ ), 1 );
        v( tt & ti ) = 4;
        tt = tt & ~ti;
    end
    v( tt ) = sign( x.basis_( 1, tt ) ) + 2;
end

tt = ~tt & ~isnan( v );
if any( tt ),
    temp = v( tt );
    temp = temp + 9;
    v( tt ) = temp;
    if ~isreal( x.basis_ ),
        ti = any( imag( x.basis_ ), 1 );
        v( tt & ti ) = 14;
    end
end

tt = isnan( v );
v( tt ) = 19;

if nnz( cvx___.exp_used ),
    tt = find( ( v == 19 | v == 12 | v == 13 ) & q == 1 );
    if ~isempty( tt ),
        [ rx, cx, vx ] = find( x.basis_( :, tt ) );
        qq = reshape( cvx___.logarithm( rx ), size( vx ) ) & ( vx > 0 );
        v( tt( cx( qq ) ) ) = 16 + sign( cvx___.vexity( cvx___.logarithm( rx( qq ) ) ) );
    end
    tt = find( ( v == 12 | v == 13 ) & q > 1 );
    if ~isempty( tt ),
        [ rx, cx, vx ] = find( x.basis_( :, tt ) );
        qq = ( ~reshape( cvx___.logarithm( rx ), size( vx ) ) & ( rx > 1 ) ) | vx < 0;
        tt( cx( qq ) ) = [];
        v( tt ) = 18;
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
