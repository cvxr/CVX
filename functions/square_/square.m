function y = square( x )

%SQUARE    Square.
%   SQUARE(X) is an array of the same size as X, whose elements are the
%   squares of the elements of X.
%
%   Disciplined convex programming information:
%       If X is real, then SQUARE(X) is convex and nonmonotonic in X. If X
%       is complex, then SQUARE(X) is neither convex nor concave. Thus when
%       when use in CVX expressions, X must be real and affine.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' }, { 'l_valid' }, ...
        { 'r_affine', 'p_convex', 'n_concave' } );
    P.funcs = { @square_cnst, @square_logv, @square_affn };
end

try
    y = cvx_unary_op( P, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
end

function y = square_cnst( x )
y = builtin( 'power', x, 2 );

function y = square_logv( x )
y = exp( 2 * log( x ) );

function y = square_affn( x ) %#ok
cvx_begin
    epigraph variable y( size(x) ) nonnegative_
    { linearize_abs( x ), 1, y } == rotated_lorentz( size(x), 0 ); %#ok
cvx_end
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
