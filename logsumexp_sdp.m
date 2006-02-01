function cvx_optval = logsumexp_sdp( x, dim, tol )

%LOGSUMEXP_SDP    SDP-based approximation of log(sum(exp(x))).
%
%   LOGSUMEXP_SDP(X) computes an approximation of the function
%   LOG(SUM(EXP(X))) using semidefinite programming techniques.
%   The approximation is chosen so that, to within the numerical
%   tolerance of the SDP solver,
%         Y <= LOGSUMEXP_SDP(X) <= Y + TOL
%   where Y = LOG(SUM(EXP(X))).
%
%   Our specific choice of a one-sided, absolute approximation has
%   two important consequences. First of all, the one-sidedness insures
%   that constraints that utilize LOGSUMEXP_SDP are conservative; that
%   is, they are tighter than if they were to use LOGSUMEXP exactly.
%   So if the approximate disciplined convex program is feasible, so
%   is the original.
%
%   Secondly, for geometric programs, the absolute tolerance TOL
%   translates to a *relative* tolerance of EXP(TOL) in a posynomial
%   constraint. That is, given a constraint
%       P(X) <= M(X)
%   where P(X) is posynomial and M(X) is monomial, the approximation
%   has the effect of enforcing a constraint
%       P(X) <= M(X)/(1+E)
%   for some unknown 0 <= E <= TOL, if TOL is sufficiently small.
%
%   If X is a matrix, LOGSUMEXP_SDP(X) will perform its computations
%   along each column of X. If X is an N-D array, LOGSUMEXP_SDP(X)
%   will perform its computations along the first dimension of size
%   other than 1. LOGSUMEXP_SDP(X,DIM) will perform its computations
%   along dimension DIM.
%
%   LOGSUMEXP_SDP(X,[],TOL) and LOGSUMEXP_SDP(X,DIM,TOL) allow you to
%   specify a different tolerance level TOL. The function will attempt
%   to select a polynomial that guarantees that
%         Y <= LOGSUMEXP_SDP(X) <= Y + TOL.
%   where Y = LOG(SUM(EXP(X))). Note that a fixed set of polynomials
%   have been hard-coded into this function. So if TOL is too small, an
%   error will result. In particular, the tightest tolerance currently
%   available is approximately TOL = +3.71E-005 * SIZE(X,DIM).
%
%   Disciplined convex programming information:
%       LOGSUMEXP_SDP(X) is convex an nondecreasing in X; therefore, X
%       must be convex (or affine).

if ~isreal( x ),
    error( 'Input must be real.' );
end
sx = size( x );

if nargin < 2 | isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim, true ),
    error( 'Second argument must be a valid dimension.' );
end

if nargin < 3 | isempty( tol ),
    tol = 0.01;
elseif ~isnumeric( tol ) | length( tol ) ~= 1 | ~isreal( tol ) | tol <= 0 | tol >= 1,
    error( 'tol must be a number between 0 and 1, exclusuve.' );
end

if length( sx ) < dim,
    sx( end + 1 : dim ) = 1;
end
nx = sx( dim );
if nx == 0,
    sx( dim ) = 1;
    cvx_optval = -Inf * ones( sx );
elseif any( sx == 0 ),
    cvx_optval = zeros( sx );
end

persistent polynomials tolerances offsets tols_lse2 xmax_lse2 poly_lse2;
if isempty( offsets ),
    tolerances = [ ...
       +1.096784836686553e-002 ...
       +1.584131350487621e-003 ...
       +2.389330716917604e-004 ...
       +3.701856466820086e-005 ...
    ];
    offsets = [ ...
       +4.594230245949757e+000 ...
       +6.482487543577793e+000 ...
       +8.370930684437980e+000 ...
       +1.022866382941549e+001 ...
    ];
    polynomials = { ...
        [ +2.198615202645529e+000, -2.086429903560452e+000, +8.792771813671837e-001, -2.430328819132307e-003, +1.096784836686576e-002 ], ...
        [ +5.458222035373999e+000, -1.051861889768797e+001, +8.517456747499884e+000, -2.910507777096101e+000, +4.511973992112296e-001, +6.663613484718728e-004, +1.584131350488056e-003 ], ...
        [ +1.410366857538211e+001, -4.131180347250113e+001, +5.204111132410120e+001, -3.410477327047925e+001, +1.235555041712883e+001, -2.274261169013964e+000, +1.906743183323483e-001, -4.056560218423852e-004, +2.389330716914037e-004 ], ...
        [ +3.626227994911397e+001, -1.423085954538006e+002, +2.456795520503084e+002, -2.371935030334587e+002, +1.389042278469432e+002, -4.985702916188342e+001, +1.068153251732325e+001, -1.232933665173396e+000, +6.452568436649889e-002, -9.375230384581923e-005, +3.701856466438633e-005 ], ...
    };
    tols_lse2 = [ ...
       +2.604468524197358e-003 ...
    ];
    xmax_lse2 = [ ...
       +3.022262664158867e+000 ...
    ];
    poly_lse2 = { ...
        [ 1.316208971334575e+000, -4.308372321551013e+000, 5.360536986239956e+000, -3.716409635631861e-002, 6.936575930158925e-001 ], ...
    };
end

%
% Determine the computation method. Note that convex inputs cannot use the lse2()
% method, because it uses an intermediate computation that is nonmonotonic.
%

nlevs  = ceil(log2(nx));
lintol = - log( 1 - nx * tolerances );
dectol = - nlevs * log( 1 - 2 * tolerances );
ls2tol = + nlevs * tols_lse2;
degs = [ min([find(ls2tol<=tol),Inf]), min([find(dectol<=tol),Inf]), min([find(lintol<=tol),Inf]) ];
if ~isnumeric( x ) & ~cvx_isaffine( x ), 
    degs(1) = Inf; 
end
if all( isinf( degs ) ),
    if 1,
        tmax = min( [ lintol(end), dectol(end) ] );
    else,
        tmax = min( [ lintol(end), ls2tol(end), dectol(end) ] );
    end
    error( sprintf( 'A polynomial of required accuracy (%g) has not been supplied.\nConsider raising the tolerance to %g or greater to proceed.', tol, tmax ) );
end
nnx = nx;
npairs = 0;
for k = 1 : nlevs,
    npairs = npairs + floor( 0.5 * nnx );
    nnx = ceil( 0.5 * nnx );
end
cplx = ( 0.5 .* ( degs + 2 ) .* ( degs + 3 ) + [ 5, 2, 2 ] ) .* [ npairs, 2 * npairs, nx ];
[ cmin, dndx ] = min( cplx );
use_lse2 = dndx == 1;
if use_lse2,
    xoff = xmax_lse2(degs(1));
    p = poly_lse2{degs(1)};
else,
    xoff = offsets(degs(dndx));
    p = polynomials{degs(dndx)};
    if dndx == 3, npairs = 1; end
end

%
% Permute the matrix, if needed, so the geometric mean can be taken
% along the first dimension.
%

if dim > 1 & any( sx( 1 : dim - 1 ) > 1 ),
    perm = [ dim, 1 : dim - 1, dim + 1 : length( sx ) ];
    x = permute( x, perm );
    sx = sx( perm );
    dim = 1;
else,
    perm = [];
end
nv = prod( sx ) / nx;
x = reshape( x, nx, nv );

%
% Perform the computation.
%

cvx_begin sdp separable
    variable y( 1, nv )
    minimize( y )
    if npairs > 1,
        variable xtemp( npairs - 1, nv );
        xq = reshape( [ x ; xtemp ], 2, npairs * nv );
        yq = reshape( [ xtemp ; y ], 1, npairs * nv );
    else,
        xq = x;
        yq = y;
    end
    if use_lse2,
        variables w( 1, npairs * nv ) v( 1, npairs * nv )
        abs( [0.5,-0.5]*xq ) <= w + v;
        w <= xoff;
        w >= 0;
        v >= 0;
        polyval_sdp( p, w / xoff ) + v + [0.5,0.5]*xq <= yq;
    else,
        xy = xq - ones(size(xq,1),1) * yq;
        xy = cvx_accept_convex( max( 0, xy + xoff ) / xoff );
        sum( polyval_sdp( p, xy ), 1 ) <= 1;
    end
cvx_end

%
% Reverse the reshaping and permutation steps
%

sx( dim ) = 1;
cvx_optval = reshape( cvx_optval, sx );
if ~isempty( perm ),
    cvx_optval = ipermute( cvx_optval, perm );
end
