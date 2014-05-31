function y = pdom( x )

% PDOM    Positive domain.
%    PDOM(X)=X, with a domain constraint that X >= 0. X must be real.
%
%    Like POS, PDOM can be useful as an argument to functions that require
%    a positive argument. For instance, for affine X, PDOM(X)^3 is X^3 when
%    X is positive, and constrains X to be positive as well. In contrast,
%    POS(X)^3=MAX(X,0)^3 is also X^3 when X is positive; but it does not
%    constrain X, and instead returns 0 when X is negative. 
%
%    One useful application of PDOM(X) is to implement 1/X for positive X:
%    that is, 1/PDOM(X) = INV_POS(X) = POW_P(X,-1).
%
%    Disciplined convex programming information:
%        From a strict mathematical standpoint, PDOM(X) is best interpreted
%        as a concave increasing function of X. A simple, sign-independent
%        application of the DCP ruleset, therefore, requires that X be 
%        concave or real affine; convex arguments are rejected. Under
%        sign-sensitive analysis, however, PDOM(X) is best understood as
%        an *identify* function PDOM(X)=X that also adds the constraint
%        that X >= 0. So, for instance:
%        --- If X is nonnegative, then PDOM(X)=X, and even convex arguments
%            are accepted. So, for instance, PDOM(SQUARE(X))=SQUARE(X).
%        --- If X is convex with unknown or negative sign, then PDOM(X) 
%            returns a DCP error, as expected.
%        --- Otherwise, PDOM(X)=X, and a constraint X>=0 is added to the 
%            model. The curvature of PDOM(X) is the same as X: that is, if
%            X is affine, then PDOM(X) is nonnegative affine; if X is
%            concave, then PDOM(X) is nonnegative concave.
%        If X is nonpositive, then a warning will be issued, since the only
%        potentially feasible value of X will be X == 0.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'negative' }, { 'nonnegative' }, ...
        { 'p_nonconst' }, { 'n_concave' }, { 'r_affine' }, ...
        { 'concave' }, [0,1,2,3,4,5] );
    P.funcs = { @pdom_nneg, @pdom_nneg, @pdom_npos, @pdom_affn, @pdom_cncv };
end
y = cvx_unary_op( P, x );

function y = pdom_nneg( x )
% Nonnegative
y = x;

function y = pdom_npos( x )
% Nonpositive
warning( 'CVX:Warning', ...
    [ 'Disciplined convex programming warning:\n', ...
      'Almost certainly infeasible: pdom( {%s} ).' ], ...
    cvx_class( x, true, true, true ) );
y = cvx( size(x), [] );
cvx_begin
    x >= 0; %#ok
cvx_end

function y = pdom_affn( x ) %#ok
cvx_begin set
    variable y(size(x)) nonnegative
    y == x; %#ok
cvx_end

function y = pdom_cncv( x ) %#ok
cvx_begin set
    hypograph variable y(size(x)) nonnegative
    y <= x; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
