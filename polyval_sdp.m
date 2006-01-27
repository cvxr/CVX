function cvx_optval = polyval_sdp( p, x )

%POLYVAL_SDP    SDP-based evaluation of polynomials.
%
%   POLYVAL_SDP uses a semidefinite program to compute the value of
%   the polynomial represented by the vector p. The format of p is
%   identical to that required by the standard POLYVAL function.
%   The vector x must be real.
%
%   If p is convex and bounded below, then POLYVAL_SDP( p, x )
%   produces the same result as POLYVAL( p, x ) to within the
%   numerical tolerance of the CVX solver. Otherwise, it returns
%   the value of the function whose epigraph is the convex hull of
%   the polynomial. In particular, if p is unbouned below, then
%   POLYVAL_SDP( p, x ) will return -Inf. Thus for meaningful
%   results, the degree of the polynomial must be even; i.e.,
%   LENGTH( p ) must be odd.
%
%   If x is an array, the result will be computed elementwise.
%
%   Disciplined convex programming information:
%       POLYVAL_SDP(P,X) is convex and nonmotonic in X; therefore, in CVX
%       expressions, X must be affine. P must always be constant.

%
% Check the polynomial
%

n = length( p );
if isempty( p ) | ~isa( p, 'double' ) | ~isreal( p ) | numel( p ) ~= n | any( isinf( p ) | isnan( p ) ),
    error( 'First argument must be a non-empty real vector.' );
end

%
% Check the input
%

sx = size( x );
if ~isreal( x ),
    error( 'Second argument must be real.' );
elseif ~cvx_isaffine( x ),
    error( sprintf( 'Disciplined convex programming error:\n   POLYVAL_SDP(P,X) requires that X be affine.' ) );
end

%
% Quick exits
%

n = min( find( p ) );
if isempty( n ),
    cvx_optval = zeros( sx );
    return
end
p( 1 : n - 1 ) = [];
n = length( p );
if rem( n, 2 ) ~= 1 | p( 1 ) < 0,
    cvx_optval = -Inf*ones( sx );
    return
end

%
% Build and solve the SDP
%

degr = n - 1;
deg2 = 0.5 * degr + 1;
nv   = prod( sx );
p    = p(:).';
cvx_begin sdp separable
    variable y(sx);
    variable P(deg2,deg2,sx) hankel;
    P >= 0;
    1 == P(1,1,:);
    x == reshape( P(2,1,:), sx );
    y == reshape( p(end:-1:1) * [ reshape( P(1,:,:), deg2, nv ) ; reshape( P(2:end,end,:), deg2-1, nv ) ], sx );
    minimize( y );
cvx_end

