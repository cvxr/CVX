function y = pdom( x )

% PDOM(X) internal CVX implementation.

persistent remap funcs
if isempty( funcs ),
    remap = max( 0, cvx_remap( ...
        { 'negative' }, ...
        { 'nonnegative', 'p_nonconst' }, ...
        { 'n_concave' }, ...
        { 'r_affine' }, ...
        { 'concave' }, [-1,1,2,3,4] );
    funcs = { @pdom_1, @pdom_2, @pdom_3, @pdom_4 };
end

try
    y = unary_op( 'pdom', funcs, remap, x );
catch exc
    throw( exc );
end

function y = pdom_1( x )
% Nonnegative
y = x;

function y = pdom_2( x )
% Nonpositive
sx = x.size_;
y = cvx( sx, [] );
y( 'Disciplined convex programming warning:\n    Almost certainly infeasible: p( {%s} ).', ...
    cvx_class( x, true, true, true ) );
cvx_begin
    x >= 0; %#ok
cvx_end

function y = pdom_3( x )
% Real affine
y = [];
sx = x.size_;
cvx_begin
    variable y( sx ) nonnegative
    y == x; %#ok
cvx_end

function y = pdom_4 ( x )
% Concave
y = [];
sx = x.size_;
cvx_begin
    hypograph variable y( sx ) nonnegative
    y <= x; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.