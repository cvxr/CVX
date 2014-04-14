function y = abs( x )

%Disciplined convex/geometric programming information for ABS:
%   ABS(X) is convex and nonmonotonic in X. Therefore, according to
%   the strict rules of DCP, X must be affine. However, because of
%   its special structure, CVX considers the sign of X as well. So,
%   for instance, if X is known to be nonnegative, then ABS(X)=X.

persistent remap funcs
if isempty( remap ),
    remap = cvx_remap( ...
        { 'constant' }, ...
        { 'p_nonconst' }, ...
        { 'n_nonconst' }, ...
        { 'r_affine' }, ...
        { 'c_affine' } );
    funcs = { @abs_1, @abs_2, @abs_3, @abs_4, @abs_5 };
end

try
    y = unary_op( 'abs', funcs, remap, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = abs_1( x )
% Constant
y = cvx( abs( cvx_constant( x ) ) );

function y = abs_2( x )
% Positive any
y = x;

function y = abs_3( x )
% Negative any
y = -x;

function y = abs_4( x )
% Real affine
y = [];
sx = x.size_;
cvx_begin
    epigraph variable y( sx )
    { x, y } == lorentz( sx, 0 ); %#ok
    cvx_setnneg( y );
cvx_end

function y = abs_5( x )
% Complex affine
y = [];
sx = x.size_;
cvx_begin
    epigraph variable yt( sx )
    { x, y } == complex_lorentz( sx, 0 ); %#ok
    cvx_setnneg( y );
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
