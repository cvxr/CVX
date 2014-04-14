function y = linearize( x )

% LINEARIZE    Linearize.
%    For real affine X, Y = X. 
%    For convex X, Y is linear, and satisfies Y >= X.
%    For concave X, Y is linear, and satisfies Y <= X.
% This is used primarily within CVX functions to efficiently implement
% certain monotonic functions.

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { 'real' }, ...
        { 'r_affine' }, ...
        { 'convex' }, ...
        { 'concave' } );
    params.constant = 1;
    params.funcs = { @linearize_1, @linearize_1, @linearize_2, @linearize_3 };
    params.name = 'linearize';
end

try
    y = cvx_unary_op( params, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = linearize_1( x )
y = x;

function y = linearize_2( x ) %#ok
cvx_begin
    variable y(size(x))
    x <= y; %#ok
cvx_end

function y = linearize_3( x )
cvx_begin
    variable y(size(x))
    x >= y; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.