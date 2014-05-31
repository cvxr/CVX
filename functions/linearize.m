function y = linearize( x )

% LINEARIZE    Linearize.
%    For real affine X, Y = X. 
%    For convex X, Y is linear, and satisfies Y >= X.
%    For concave X, Y is linear, and satisfies Y <= X.
% This is used primarily within CVX functions to efficiently implement
% certain monotonic functions.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'r_affine' }, ...
        { 'convex' }, { 'concave' }, [2,3,4] );
    P.funcs = { [], @lin_affn, @lin_cnvx, @lin_cncv };
end
y = cvx_unary_op( P, x );

function y = lin_affn( x )
y = x;

function y = lin_cnvx( x ) %#ok
cvx_begin set
    variable y(size(x))
    x <= y; %#ok
cvx_end

function y = lin_cncv( x ) %#ok
cvx_begin set
    variable y(size(x))
    x >= y; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
