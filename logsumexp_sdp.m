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
%   has the effect of, in the worst case, producing
%       P(X) <= M(X)*EXP(-TOL).
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
if sx( dim ) == 0,
    sx( dim ) = 1;
    cvx_optval = -Inf * ones( sx );
elseif any( sx == 0 ),
    cvx_optval = zeros( sx );
end

ntol = ( 1 - exp( -tol ) ) / sx( dim );
if ntol >= +1.096784836686553e-002,
    xoff = +4.594230245949757e+000;
    p = [ +2.198615202645529e+000 ...
            -2.086429903560452e+000 ...
            +8.792771813671837e-001 ...
            -2.430328819132307e-003 ...
            +1.096784836686576e-002 ];
elseif ntol >= 1.584131350487621e-003,
    xoff = +6.482487543577793e+000;
    p = [ +5.458222035373999e+000 ...
            -1.051861889768797e+001 ...
            +8.517456747499884e+000 ...
            -2.910507777096101e+000 ...
            +4.511973992112296e-001 ...
            +6.663613484718728e-004 ...
            +1.584131350488056e-003 ];
elseif ntol >= +2.389330716917604e-004,
    xoff = +8.370930684437980e+000;
    p = [ +1.410366857538211e+001 ...
            -4.131180347250113e+001 ...
            +5.204111132410120e+001 ...
            -3.410477327047925e+001 ...
            +1.235555041712883e+001 ...
            -2.274261169013964e+000 ...
            +1.906743183323483e-001 ...
            -4.056560218423852e-004 ...
            +2.389330716914037e-004 ];
elseif ntol >= +3.701856466820086e-005,
    xoff = +1.022866382941549e+001;
    p = [ +3.626227994911397e+001 ...
            -1.423085954538006e+002 ...
            +2.456795520503084e+002 ...
            -2.371935030334587e+002 ...
            +1.389042278469432e+002 ...
            -4.985702916188342e+001 ...
            +1.068153251732325e+001 ...
            -1.232933665173396e+000 ...
            +6.452568436649889e-002 ...
            -9.375230384581923e-005 ...
            +3.701856466438633e-005 ];
else,
    tmax = -log( 1 - 3.701856466820086e-005 * sx(dim) );
    error( sprintf( 'A polynomial of required accuracy (%g) has not been supplied.\nConsider reducing the tolerance to %g to proceed.', tol, tmax ) )';
end

sy = sx;
sy( dim ) = 1;
ndxs = cell( 1, length( sx ) );
[ ndxs{:} ] = deal( ':' );
ndxs{dim} = ones( 1, sx( dim ) );

cvx_begin sdp
    variable y(sy);
    if cvx_isconstant( x ),
        minimize sum( y( : ) );
    else,
        minimize y;
    end
    temp = ( ( x - y( ndxs{:} ) ) + xoff ) / xoff;
    temp = max( temp, 0 );
    temp = cvx_accept_convex( temp );
    sum( polyval_sdp( p, temp ), dim ) <= 1;
cvx_end
if isnumeric( x ),
    cvx_optval = y;
end

