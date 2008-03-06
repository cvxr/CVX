function y = pow_abs( x, p )

%POW_ABS   Internal cvx version.

error( nargchk( 2, 2, nargin ) );
if ~cvx_isconstant( p ),
    error( 'Second argument must be constant.' );
elseif ~isreal( p ),
    error( 'Second argument must be real.' );
end
p = cvx_constant( p );
if nnz( isinf( p ) | isnan( p ) ),
    error( 'Inf or NaN not permitted here.' );
elseif nnz( p < 1 ),
    error( 'Second argument must be greater than or equal to 1.' );
end

%
% Check sizes
%

sx = size( x ); xs = all( sx == 1 );
sp = size( p ); ps = all( sp == 1 );
if xs,
    sy = sp;
elseif ps | isequal( sx, sp ),
    sy = sx;
else
    error( 'Matrix dimensions must agree.' );
end

%
% Determine the expression types
% 0 : invalid
% 1 : constant .^ p
% 2 : real affine, complex affine .^ 1
% 3 : real affine .^ (p>1)
% 4 : complex affine .^ (p>1)
%

persistent remap
if isempty( remap ),
    remap_x1 = cvx_remap( 'constant' );
    remap_x2 = cvx_remap( 'real-affine' );
    remap_x3 = cvx_remap( 'affine' );
    remap_x4 = cvx_remap( 'log-valid' ) & ~remap_x1;
    remap_x3 = remap_x3 & ~remap_x2;
    remap_x2 = remap_x2 & ~remap_x1;
    remap    = [1;1] * remap_x1 + [2;6] * remap_x2 + [2;3] * remap_x3 + [4;5] * remap_x4;
end
v = 1 + ( p > 1 );
v = remap( v(:)' + 2 * ( cvx_classify( x ) - 1 ) );
t = v == 6;
if any( t ) & ~ps,
    [ pk, pi, pj ] = unique( p( t ) );
    v( t ) = pj + 5;
end

%
% Process each type of expression one piece at a times
%

xt = x;
pt = p;
vu = unique( v );
nv = length( vu );
if nv ~= 1,
    if cvx_isconstant( x ),
        y = zeros( sy );
    else
        y = cvx( sy, [] );
    end
end
for k = 1 : nv,
    
    %
    % Select the category of expression to compute
    %
    
    vk = vu( k );
    if nv ~= 1,
        t = v == vk;
        if ~xs, xt = cvx_subsref( x, t ); end
        if ~ps, pt = cvx_subsref( p, min(find(t)) ); end
    end
    sz = size(xt);
    nd = length(sz)+1;
    
    %
    % Perform the computations
    %
    
    switch vk,
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: pow_pos( {%s}, p )', cvx_class( xt, true, true ) ) ); 
        case 1,
            % Constant .^ p
            yt = abs( cvx_constant( xt ) ) .^ pt;
        case 2,
            % affine .^ 1
            yt = abs( xt );
        case 3,
            % Complex affine .^ p
            cvx_begin
                epigraph variable yt(sz)
                geomean( cat( nd, yt, ones(sz) ), nd, cvx_geomean_map( pt ), true, 'cabs' ) == xt;
            cvx_end
        case 4,
            % Log-valid .^ 1
            yt = xt;
        case 5,
            % Log-valid .^ p
            yt = xt .^ pt;
        otherwise,
            % Real affine .^ p
            cvx_begin
                epigraph variable yt(sz)
                geomean( cat( nd, yt, ones(sz) ), nd, cvx_geomean_map( pt ), true, 'abs' ) == xt;
            cvx_end
    end
    
    %
    % Store the results
    %
    
    if nv ~= 1,
        y = cvx_subsasgn( y, t, yt );
    elseif isequal(sz,sy),
        y = yt;
    else
        y = yt * ones(sy);
    end
    
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
