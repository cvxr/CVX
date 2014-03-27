function v = cvx_classify( x )

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - nonpositive concave
% 6  - concave
% 7  - real affine
% 8  - convex
% 9  - nonnegative convex
% 10 - complex affine
% 11 - log concave
% 12 - log affine
% 13 - log convex monomial
% 14 - log convex posynomial
% 15 - invalid

global cvx___
v = full( cvx_vexity( x ) );
v = reshape( v, 1, prod( x.size_ ) );
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
    temp = temp + 7;
    v( tt ) = temp;
    if ~isreal( x.basis_ ),
        ti = any( imag( x.basis_ ), 1 );
        v( tt & ti ) = 10;
    end
end

tt = isnan( v );
v( tt ) = 15;

if nnz( cvx___.exp_used ),
    tt = find( ( v == 15 | v == 8 | v == 9 ) & q == 1 );
    if ~isempty( tt ),
        [ rx, cx, vx ] = find( x.basis_( :, tt ) );
        qq = reshape( cvx___.logarithm( rx ), size( vx ) ) & ( vx > 0 );
        v( tt( cx( qq ) ) ) = 12 + sign( cvx___.vexity( cvx___.logarithm( rx( qq ) ) ) );
    end
    tt = find( ( v == 8 | v == 9 ) & q > 1 );
    if ~isempty( tt ),
        [ rx, cx, vx ] = find( x.basis_( :, tt ) );
        qq = ( ~reshape( cvx___.logarithm( rx ), size( vx ) ) & ( rx > 1 ) ) | vx < 0;
        tt( cx( qq ) ) = [];
        v( tt ) = 14;
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
