function y = square( x )

% SQUARE   Internal cvx version.

persistent remap funcs
if isempty( funcs ),
    remap = cvx_remap( ...
        { 'real' }, ...
        { 'l_valid' }, ...
        { 'r_affine' }, ...
        { 'p_convex', 'n_concave' } );
    funcs = { @square_1, @square_2, @square_3, @square_4 };
end

try
    y = unary_op( 'square', funcs, remap, x );
catch exc
    throw( exc );
end

function y = square_1( x )
y = cvx( cvx_constant( x ) .^ 2 );

function y = square_2( x )
y = exp( 2 * log( x ) );

function y = square_3( x )
y = [];
sx = x.size_;
cvx_begin
    epigraph variable y( sx )
    { x, 1, y } == rotated_lorentz( sx, length(sx)+1 ); %#ok
    cvx_setnneg(y);
cvx_end

function y = square_4( x )
y = [];
sx = x.size_;
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
