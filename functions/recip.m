function y = recip( x )

% RECIP RECIP(X) = 1.0 / X, except X must be real and nonzero.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'positive', 'negative' }, { 'l_valid' }, ...
        { 'p_concave' }, { 'n_convex' } );
    P.funcs = { @recip_cnst, @recip_logv, @recip_posc, @recip_negc };
end

try
    y = cvx_unary_op( P, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = recip_cnst( x )
y = 1.0 ./ x;

function y = recip_logv( x )
y = exp( -log( x ) );

function y = recip_posc( x ) %#ok
% Positive concave
sx = size( x );
cvx_begin
    epigraph variable y( sx ) nonnegative_
    { 1, linearize(x), y } == rotated_lorentz( sx, 0 ); %#ok
cvx_end

function y = recip_negc( x ) %#ok
% Negative convex    
y = -recip_posc( -x );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
