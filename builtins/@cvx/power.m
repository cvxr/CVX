function z = power( x, y, op )

%POWER   Internal cvx version.

persistent bparam
if isempty( bparam ),
    bparam.map = cvx_remap( ...
        { { 'constant' } }, ...
        { { 'positive' }, { 'convex', 'concave' } }, ...
        { { 'valid' }, { 'real' } }, [1,2,3] );
    bparam.funcs = { @power_1, @power_2, @power_3 };
end

try
    bparam.name = op;
    z = cvx_binary_op( bparam, x, y );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function z = @power_1( x, y )
z = cvx( power( cvx_constant( x ), cvx_constant( y ) ) );

function z = @power_2( x, y )
z = exp( log( cvx_constant( x ) ) .* y );

function z = @power_3( x, y )
nx = x.size_(1);
ny = numel( y );
if ny == 1
    pu = y;
    multp = false;
else
    pu = sort( y' ); %#ok'
    pu = pu( [ true, diff(pu) ~= 0 ] );
    multp = length( pu ) ~= 1;
end
if multp,
    zp = cvx( sz, [] );
    vp = vx;
    xp = x;
    yp = y;
end
for pk = pu,
    if multp,
        tp = yp == pk;
        if nx > 1, 
            x = xp.basis_( :, tp );
            x = cvx( [size(x,2),1], x );
            vx = vp( tp ); 
        end
        if ny > 1, 
            y = yp( tp ); 
        end
    end
    np = numel(y);
    y1 = y(1);
    if y1 > 1,
        z = power_cvx( x, vx, y1 );
    elseif y1 == 1,
        z = x;
    elseif y1 > 0,
        z = power_ccv( x, vx, y1 );
    elseif y1 == 0,
        z = cvx( ones(size(x)) );
    else
        z = power_neg( x, vx, y1 );
    end
    if islogical( z ),
        errs(end+1,:) = { cvx_subsref( x, z ), y }; %#ok
        continue
    end
    if nx == 1 && np > 1,
        z.size_(1) = np;
        z.basis_ = repmat( z.basis_, [1,np] );
    end
    if multp,
        zp.basis_(1:size(z.basis_,1),tp) = z.basis_;
    end
end
if multp
    z = zp;
end
if ~isempty( errs ),
    throw( cvx_dcp_error( errs, op ) );
elseif mult
    z = zm;
else
    z.size_ = sz;
end

%
% P > 1, integer: X affine, p-convex, n-concave
% P > 1, non-integer: X affine, p-convex
%

function y = power_cvx( x, v, p )
persistent remap remap_i remap_e
if isempty( remap_e ),
    remap   = cvx_remap( 'affine', 'p_convex' );
    remap_i = cvx_remap( 'p_convex', 'n_concave' );
    remap_e = remap_i | remap;
end
if ~isempty( v ),
    isint = rem( p, 1 ) == 0;
    if isint,
        isevn = rem( p, 2 ) == 0;
        if isevn,
            v = remap_e( v );
        else
            v = remap_i( v );
        end
    else
        v = remap( v );
    end
    if ~all( v ),
        y = v == 0;
        return
    end
    v = 1 - 2 * ( v == 3 );
    x = cvx_accept_convex( v .* x );
else
    isint = false;
end
ne = 0;
nx = numel( x );
while rem(p,2) == 0,
    ne = ne + 1;
    p = p * 0.5;
end
cvx_begin
    epigraph variable y(nx)
    if ne > 0,
        variables z2(nx,ne-1)
        if p == 1,
            w = y; %#ok
        else
            variable w(nx)
        end
        { [x,z2], 1, [z2,w] } == rotated_lorentz( [ nx, ne ], 3 ); %#ok
    else
        w = x;
    end
    if p > 1,
        { [y,ones(nx,1)], w } == geo_mean_cone( [ nx, 2 ], 2, [1/p,1-1/p], 'func' );  %#ok
    end
    cvx_setnneg(y)
cvx_end
if isint,
    y = v .* y; 
end

%
% 0 < P < 1: X concave
%

function y = power_ccv( x, v, p )
persistent remap
if isempty( remap ),
    remap = ~cvx_remap( 'concave' );
end
v = remap( v );
if any( v ), 
    y = v;
    return
end
p = 1.0 / p;
if p ~= floor( p ),
    [ nn, dd ] = rat( p );
    p = nn / dd;
end    
cvx_begin
    hypograph variable y(numel(x))
    power_cvx( y, [], p ) <= x; %#ok
    cvx_setnneg(y);
cvx_end

%
% P < 0, integer: p-concave, n-convex
% P < 0, non-integer: p-concave
%

function y = power_neg( x, v, p )
persistent remap remap_i
if isempty( remap_i ),
    remap   = cvx_remap( 'p_concave' );
    remap_i = cvx_remap( { 'p_concave' }, { 'n_convex' }, [1,-1] );
end
if rem( p, 1 ) == 0,
    % Integer: positive and negative lobes
    v = remap_i( v );
else
    % Non-integer: positive lobe only
    v = remap( v );
end
if ~all( v ),
    y = v == 0;
    return
end
nx = numel(x);
z  = cvx_accept_concave( v .* x );
cvx_begin
    epigraph variable y(nx)
    if p == 1,
        { 1, z, y } == rotated_lorentz( [nx,1], 2, 0 ); %#ok
    else
        { [z,y], 1 } == geo_mean_cone( [nx,2], 2, [-p,1], 'func' ); %#ok
    end
    cvx_setnneg(y);
cvx_end
y = v .* y;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

