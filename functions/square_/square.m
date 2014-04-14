function y = square( x )

%SQUARE    Square.
%   SQUARE(X) is an array of the same size as X, whose elements are the
%   squares of the elements of X.
%
%   Disciplined convex programming information:
%       If X is real, then SQUARE(X) is convex and nonmonotonic in X. If X
%       is complex, then SQUARE(X) is neither convex nor concave. Thus when
%       when use in CVX expressions, X must be real and affine.

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { 'real' }, ...
        { 'l_valid' }, ...
        { 'r_affine' }, ...
        { 'p_convex', 'n_concave' } );
    params.funcs = { @square_1, @square_2, @square_3, @square_4 };
    params.constant = 1;
    params.name = 'square';
end

try
    y = cvx_unary_op( params, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
end

function y = square_1( x )
y = x .^ 2;

function y = square_2( x )
y = exp( 2 * log( x ) );

function y = square_3( x ) %#ok
sx = size( x );
cvx_begin
    epigraph variable y( sx )
    { x, 1, y } == rotated_lorentz( sx, length(sx)+1 ); %#ok
    cvx_setnneg(y);
cvx_end

function y = square_4( x ) %#ok
sx = size( x );
cvx_begin
    epigraph variable y( sx )
    variable z( sx )
    { z, 1, y } == rotated_lorentz( sx, length(sx)+1 ); %#ok
    abs( x ) <= z; %#ok
    cvx_setnneg(y);
cvx_end
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
