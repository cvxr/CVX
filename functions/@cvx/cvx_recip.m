function y = cvx_recip( x )

%CVX_RECIP   Internal cvx version.

persistent remap funcs
if isempty( funcs ),
    remap = cvx_remap( ...
        { 'nonzero' }, ...
        { 'l_valid' }, ...
        { 'p_concave' }, ...
        { 'n_convex' } );
    funcs = { @recip_1, @recip_2, @recip_3, @recip_4 };
end

y = unary_op( '1.0 ./', funcs, remap, x );


function y = recip_1( x )
% Non-zero constant
y = cvx( 1.0 ./ cvx_constant( x ) );

function y = recip_2( x )
% Monomial, posynomial
y = exp( -log( x ) );

function y = recip_3( x )
% Positive concave
y = [];
sx = x.size_;
cvx_begin
    epigraph variable y( sx )
    variable z( sx )
    { 1, z, y } == rotated_lorentz( sx, length(sx)+1, 0 ); %#ok
    z <= x; %#ok
    cvx_setnneg(y);
cvx_end

function y = recip_4( x )
% Negative convex
y = [];
sx = x.size_;
cvx_begin
    hypograph variable y( sx )
    variable z( sx )
    { 1, z, -y } == rotated_lorentz( sx, length(sx)+1, 0 ); %#ok
    z <= -x; %#ok
    cvx_setnpos(y);
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.tx for full copyright information.
% The command 'cvx_where' will show where this file is located.
