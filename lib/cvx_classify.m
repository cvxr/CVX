function v = cvx_classify( x )

% Classifications:
% 1  - negative constant
% 2  - zero
% 3  - positive constant
% 4  - complex constant
% 5  - concave
% 6  - real affine
% 7  - convex
% 8  - complex affine
% 9  - log concave
% 10 - log affine
% 11 - log convex monomial
% 12 - log convex posynomial
% 13 - invalid

n = prod( size( x ) );
x = reshape( full( x ), 1, n );
if isreal( x ),
    v = sign( x ) + 2;
else,
    t = imag( x ) == 0;
    v = zeros( 1, n );
    v( t ) = sign( x( t ) ) + 2;
    v( ~t ) = 4;
end
v( isinf( x ) | isnan( x ) ) = 13;


