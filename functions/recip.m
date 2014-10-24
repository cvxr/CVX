function y = recip( x )

% RECIP RECIP(X) = 1.0 / X, except X must be real and nonzero.

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        { 'g_posynomial', 'gg_monomial' }, ...
        { 'positive', 'negative' }, { 'l_valid' }, ...
        { 'p_concave' }, { 'n_convex' }, [0,1,2,3,4] );
    P.funcs = { @recip_cnst, @recip_logv, @recip_posc, @recip_negc };
    P.name = '1 /';
end
y = cvx_unary_op( P, x );

function y = recip_cnst( x )
y = 1.0 ./ x;

function y = recip_logv( x )
y = exp( -log( x ) );

function y = recip_posc( x ) %#ok
% Positive concave
sx = size( x );
cvx_begin
    epigraph variable y( sx ) nonnegative_
    { sqrt(2), cvx_linearize(x), y } == rotated_lorentz( sx, 0 ); %#ok
cvx_end

function y = recip_negc( x )
% Negative convex   
y = -recip_posc( -x );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
