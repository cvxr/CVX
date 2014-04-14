function y = recip( x )

% RECIP RECIP(X) = 1.0 / X, except X must be real and nonzero.

persistent params
if isempty( params ),
    params.map = cvx_remap( ...
        { 'positive', 'negative' }, ...
        { 'l_valid' }, ...
        { 'p_concave' }, ...
        { 'n_convex' } );
    params.constant = 1;
    params.funcs = { @recip_1, @recip_2, @recip_3, @recip_4 };
    params.name = 'recip';
end

try
    y = cvx_unary_op( params, x );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = recip_1( x )
y = 1.0 ./ x;

function y = recip_2( x )
y = exp( -log( x ) );

function y = recip_3( x ) %#ok
% Positive concave
sx = x.size_;
cvx_begin
    epigraph variable y( sx )
    variable z( sx )
    { 1, z, y } == rotated_lorentz( sx, length(sx)+1, 0 ); %#ok
    z <= x; %#ok
    cvx_setnneg(y);
cvx_end

function y = recip_4( x ) %#ok
% Negative convex
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
