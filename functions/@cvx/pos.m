function y = pos( x )

% POS(X) internal CVX implementation.

persistent remap funcs
if isempty( funcs ),
    remap = cvx_remap( ...
        { 'real' }, ...
        { 'p_nonconst' }, ...
        { 'n_nonconst' }, ...
        { 'r_affine', 'convex' } );
    funcs = { @pos_1, @pos_2, @pos_3, @pos_4 };
end

try
    y = unary_op( 'pos', funcs, remap, x );
catch exc
    throw( exc );
end

function y = pos_1( x )
y = cvx( pos( cvx_constant( x ) ) );

function y = pos_2( x )
y = x;

function y = pos_3( x )
y = cvx( zeros( size( x ) ) );

function y = pos_4( x )
y = [];
sx = x.size_; %#ok
cvx_begin
    epigraph variable y( sx )
    x <= y; %#ok
    0 <= y; %#ok
    cvx_setnneg( y );
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
