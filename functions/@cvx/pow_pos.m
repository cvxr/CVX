function y = pow_pos( x, p )

%POW_POS   Internal cvx version.

%
% Check second argument
%

error( nargchk( 2, 2, nargin ) );
if ~cvx_isconstant( p ),
    error( 'Second argument must be constant.' );
elseif ~isreal( p ),
    error( 'Second argument must be real.' );
end
p = cvx_constant( p );
if nnz( isinf( p ) | isnan( p ) ),
    error( 'Inf or NaN not permitted here.' );
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
% 1 : nonpositive .^ (p<0)
% 2 : nonpositive .^ (p>0)
% 3 : positive or log-valid .^ anything 
% 4 : convex .^ (p>1), concave .^ (p<1)
%

persistent remap
if isempty( remap ),
    remap_x1 = cvx_remap( 'negative' );
    remap_x2 = cvx_remap( 'zero' );
    remap_x3 = cvx_remap( 'log-valid' );
    temp     = ~( remap_x1 | remap_x2 | remap_x3 );
    remap_x4 = cvx_remap( 'convex' ) & temp;
    remap_x5 = cvx_remap( 'concave' ) & temp;
    remap = [1;2;3] * remap_x1 + [1;3;3] * remap_x2 + ...
            [4;4;4] * remap_x3 + [0;5;5] * remap_x4 + [5;0;0] * remap_x5;
end
v = 1 + ( p >= 0 ) + ( p >= 1 );
v = remap( v(:)' + 3 * ( cvx_classify( x ) - 1 ) );
t = v == 5;
if any( t ) & ~ps,
    [ pk, pi, pj ] = unique( p( t ) );
    v( t ) = pj + 4;
end

%
% Process each type of expression one piece at a times
%

xt = x;
pt = p;
sz = sy;
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
        if ~xs, xt = cvx_subsref( x, t ); sz = size( xt ); end
        if ~ps, pt = cvx_subsref( p, min(find(t)) ); end
    end
    
    %
    % Perform the computations
    %
    
    switch vk,
        case 0,
            % Invalid
            error( sprintf( 'Disciplined convex programming error:\n    Illegal operation: pow_pos( {%s}, p )', cvx_class( xt, true, true ), pclass ) ); 
        case 1,
            % Nonpositive constant .^ (p<0)
            yt = +Inf * ones( sz );
        case 2,
            % Negative constant .^ (0<p<1)
            yt = -Inf * ones( sz );
        case 3,
            % Zero result
            yt = zeros( sz );
        case 4,
            % Positive constant or log-valid .^ (anything)
            yt = xt .^ pt;
        otherwise,
            if pt < 0,
                % Concave .^ (p<0)
                nd = length( sz ) + 1;
                if xs & any( sz > 1 ), 
                    xs = xs * ones(sz); 
                end
                cvx_begin
                    epigraph variable yt( sz )
                    geomean( cat( nd, xt, yt ), nd, cvx_geomean_map( p ), true ) >= 1;
                cvx_end
            elseif pt == 0,
                % Concave .^ (p==0)
                cvx_begin
                    hypograph variable yt( sz )
                    yt == 1;
                    xt >= 0;
                cvx_end
            elseif pt < 1,
                % Concave .^ (0<=p<=1)
                sxt = size(xt);
                nd = length( sxt ) + 1;
                cvx_begin
                    hypograph variable yt( sz )
                    geomean( cat( nd, xt, ones(sxt) ), nd, cvx_geomean_map( p ), true ) >= yt;
                cvx_end
            elseif pt <= 1,
                % Convex .^ (p==1)
                yt = max( xt, 0 );
            else
                % Convex .^ (p>1)
                nd = length( sz ) + 1;
                cvx_begin
                    epigraph variable yt( sz )
                    geomean( cat( nd, yt, ones(sz) ), nd, cvx_geomean_map( p ), true ) >= xt;
                cvx_end
            end
    end
    
    %
    % Store the results
    %
    
    if nv == 1,
        y = yt;
    else
        y = cvx_subsasgn( y, t, yt );
    end
    
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

