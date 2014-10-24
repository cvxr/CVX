function z = power( x, y, op )

%POWER   Internal cvx version.

persistent BP
if isempty( BP ),
    BP.map = cvx_remap( ...
        ... % Invalid combinations: zero ^ negative, posynomial ^ negative
        { { 'zero', 'g_posynomial', 'gg_monomial' }, { 'negative' } }, ...
        ... % Constant
        { { 'constant' } }, ...
        ... % constant ^ convex/concave
        { { 'positive' }, { 'convex', 'concave' } }, ...
        ... % geometric ^ real
        { { 'l_valid_' }, { 'real' } }, ...
        ... % other non-constant ^ real
        { { 'valid' }, { 'real' } }, [0,1,2,3,4] );
    BP.funcs = { @power_c, @power_e, @power_g, @power_p };
    BP.constant = 1;
end
if isnumeric( y ) && numel( y ) == 1,
    switch y,
        case   1, z = x;           return
        case   2, z = square( x ); return
        case 0.5, z = sqrt( x );   return 
        case  -1, z = recip( x );  return
    end
end
if nargin < 3, op = '.^'; end
BP.name = op;
z = cvx_binary_op( BP, x, y );

function z = power_c( x, y )
z = builtin( 'power', x, y );

function z = power_e( x, y )
z = exp( log( cvx_constant( x ) ) .* y );

function z = power_g( x, y )
z = exp( log( x ) .* y );

function z = power_p( x, y )
nx = numel( x );
ny = numel( y );
y  = cvx_constant( y );
vx = cvx_classify( x );
nz = max( nx, ny );
if ny == 1
    pu = y;
    multp = false;
else
    pu = sort( y' ); %#ok'
    pu = pu( [ true, diff(pu) ~= 0 ] );
    multp = length( pu ) ~= 1;
end
if multp,
    zp = cvx( nz, [] );
    vp = vx;
    xp = x;
    yp = y;
end
errs = {};
for pk = pu,
    if multp,
        tp = yp == pk;
        if nx > 1, 
            x = cvx_fastref( xp, tp );
            x = cvx( [size(x,2),1], x );
            vx = vp( tp ); 
        end
        if ny > 1, 
            y = yp( tp ); 
        end
    end
    y1 = y(1);
    if y == 2,
        z = square( x );
    elseif y1 > 1,
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
    if nx ~= nz,
        z = repmat( z, [nz,1] );
    end
    if multp,
        zp = cvx_fastasgn( zp, tp, z );
    end
end
if multp
    z = zp;
end
if ~isempty( errs ),
    cvx_dcp_error( errs, op );
elseif multp
    z = zm;
end

%
% P > 1, integer: X affine, p-convex, n-concave
% P > 1, non-integer: X affine, p-convex
%

function y = power_cvx_fast( x, p ) %#ok
ne = 0;
nx = numel( x );
while rem(p,2) == 0,
    ne = ne + 1;
    p = p * 0.5;
end
cvx_begin
    epigraph variable y(nx) nonnegative_
    if ne > 0,
        variables z2(nx,ne-1)
        if p == 1,
            w = y;
        else
            variable w(nx)
        end
        { [x,z2], 0.5, [z2,w] } == rotated_lorentz( [ nx, ne ], 3 ); %#ok
    else
        w = x;
    end
    if p > 1,
        { [y,ones(nx,1)], w } == geo_mean_cone( [ nx, 2 ], 2, [1/p,1-1/p], 'func' );  %#ok
    end
cvx_end

function y = power_cvx( x, v, p )
persistent remap remap_i remap_e
if isempty( remap_e ),
    remap   = cvx_remap( 'affine', 'p_convex' );
    remap_i = cvx_remap( 'p_convex' ) - cvx_remap( 'n_concave' );
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
    v = 1 - 2 * ( v(:) == 3 );
    x = cvx_linearize( v .* x );
else
    isint = false;
end
y = power_cvx_fast( x, p );
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
cvx_begin
    hypograph variable y(numel(x)) nonnegative_
    power_cvx_fast( y, 1.0 / p ) <= x; %#ok
cvx_end

%
% P < 0, integer: p-concave, n-convex
% P < 0, non-integer: p-concave
%

function y = power_neg( x, v, p )
persistent remap remap_i
if isempty( remap_i ),
    remap   = cvx_remap( 'p_concave' );
    remap_i = remap - cvx_remap( 'n_convex' );
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
z  = cvx_linearize( v .* x );
cvx_begin
    epigraph variable y(nx)
    if p == 1,
        { 2, z, y } == rotated_lorentz( [nx,1], 2, 0 ); %#ok
    else
        { [z,y], 1 } == geo_mean_cone( [nx,2], 2, [-p,1], 'func' ); %#ok
    end
    cvx_setnneg(y);
cvx_end
y = v .* y;

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

