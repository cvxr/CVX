function y = linearize_cabs( x )

% LINEARIZE_ABS    Linearize w/ complex absolute value.
%    For real affine X, Y = X. 
%    For complex affine X, Y is linear, and satisfies Y >= abs(X).
%    For nonneagtive convex X, Y is linear, and satisfies Y >= X.
%    For nonpositive concave X, Y is linear, and satisfies Y >= -X.
% This is used primarily within CVX functions to efficiently implement
% certain monotonic functions.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real', 'complex' }, { 'r_affine' }, ...
        { 'c_affine', 'p_convex', 'n_concave' } );
    P.funcs = { @lin_cabs_cplx, @lin_cabs_raff, @lin_cabs_nonl };
end

try
    y = cvx_unary_op( P, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = lin_cabs_cplx( x )
y = abs(x);

function y = lin_cabs_raff( x )
y = x;

function y = lin_cabs_nonl( x ) %#ok
cvx_begin
    variable y(size(x))
    abs(x) <= y; %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.