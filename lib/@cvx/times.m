function z = times( x, y )
error( nargchk( 2, 2, nargin ) );

%
% Check sizes
%

sx = size( x );
sy = size( y );
xs = all( sx == 1 );
ys = all( sy == 1 );
if xs,
    sz = sy;
elseif ys,
    sz = sx;
elseif ~isequal( sx, sy ),
    error( 'Matrix dimensions must agree.' );
else
    sz = sx;
end
nn = prod( sz );
zs = all( sz == 1 );

%
% Determine the computation methods
%

persistent remap
if isempty( remap ),
    % Constant result (constant * constant, zero * valid)
    remap_1 = cvx_remap( 'zero' )' * cvx_remap( 'valid' );
    remap_1 = remap_1 + remap_1';
    temp_1  = cvx_remap( 'constant' );
    remap_1 = remap_1 | ( temp_1' * temp_1 );

    % Constant * affine, Real * convex/concave/log-convex, positive * log-concave
    remap_2 = cvx_remap( 'nonzero', 'complex' )' * cvx_remap( 'affine' );
    remap_3 = cvx_remap( 'nonzero' )' * cvx_remap( 'convex', 'concave', 'log-convex' );
    remap_4 = cvx_remap( 'positive' )' * cvx_remap( 'log-concave' );
    remap_2 = ( remap_2 | remap_3 | remap_4 ) & ~remap_1;

    % Affine * constant, Convex/concave/log-convex * real, log-concave * positive
    remap_3 = remap_2';

    % Affine * affine
    temp_1  = ~temp_1;
    remap_4 = +( cvx_remap( 'real-affine' ) & temp_1 );
    remap_5 = +( cvx_remap( 'complex-affine' ) & temp_1 );
    remap_4 = ( remap_4' * remap_4 ) | ( remap_5' * remap_5 );

    % log-concave * log-concave, log-convex * log-convex
    remap_5 = +( cvx_remap( 'log-convex'  ) & temp_1 );
    remap_6 = +( cvx_remap( 'log-concave' ) & temp_1 );
    remap_5 = ( remap_5' * remap_5 ) | ( remap_6' * remap_6 );

    remap = remap_1 + ( 2 * remap_2 + 3 * remap_3 + 4 * remap_4 + 5 * remap_5 ) .* ~remap_1;
end
vx = cvx_classify( x );
vy = cvx_classify( y );
vr = remap( vx + size( remap, 1 ) * ( vy - 1 ) );
vu = unique( vr );
nv = length( vu );

%
% Process each computation type separately
%

x   = cvx( x );
y   = cvx( y );
xt  = x;
yt  = y;
xts = xs;
yts = ys;
if nv ~= 1,
    z = cvx( sz, [] );
end
for k = 1 : nv,

    %
    % Select the category of expression to compute
    %

    if nv ~= 1,
        t = vr == vu( k );
        if ~xs,
            xt = cvx_subsref( x, t );
            sz = size( xt );
        end
        if ~ys,
            yt = cvx_subsref( y, t );
            sz = size( yt );
        end
    end

    %
    % Apply the appropriate computation
    %

    switch vu( k ),
    case 0,

        % Invalid
        error( sprintf( 'Disciplined convex programming error:\n    Cannot perform the operation {%s}.*{%s}', cvx_class( xt, true, true ), cvx_class( yt, true, true ) ) );

    case 1,

        % constant .* constant
        cvx_optval = cvx_constant( xt ) .* cvx_constant( yt );

    case 2,

        % constant .* something
        xb = cvx_constant( xt );
        yb = yt.basis_;
        if ~xs,
            nn = prod( size( xb ) );
            if ys,
                xb = cvx_reshape( xb, [ 1, nn ] );
                if issparse( yb ) & ~issparse( xb ), 
                    xb = sparse( xb ); 
                end
            else
                n1 = 1 : nn;
                xb = sparse( n1, n1, xb( : ), nn, nn );
            end
        end
        cvx_optval = cvx( sz, yb * xb );

    case 3,

        % something .* constant
        xb = xt.basis_;
        yb = cvx_constant( yt );
        if ~ys,
            nn = prod( size( yb ) );
            if xs,
                yb = cvx_reshape( yb, [ 1, nn ] );
                if issparse( xb ) & ~issparse( yb ),
                    yb = sparse( yb );
                end
            else
                n1 = 1 : nn;
                yb = sparse( n1, n1, yb( : ), nn, nn );
            end
        end
        cvx_optval = cvx( sz, xb * yb );

    case 4,

        % affine .* affine
        nn = prod( sz );
        xA = xt.basis_; yA = yt.basis_;
        if xs & ~ys, xA = xA( :, ones( 1, nn ) ); end
        if ys & ~xs, yA = yA( :, ones( 1, nn ) ); end
        mm = max( size( xA, 1 ), size( yA, 1 ) );
        if size( xA, 1 ) < mm, xA( mm, end ) = 0; end
        if size( yA, 1 ) < mm, yA( mm, end ) = 0; end
        xB = xA( 1, : ); xA( 1, : ) = 0;
        yB = yA( 1, : ); yA( 1, : ) = 0;
        cyA   = conj( yA );
        alpha = sum( real( xA .* yA ), 1 ) ./ max( sum( cyA .* yA, 1 ), realmin );
        adiag = sparse( 1 : nn, 1 : nn, alpha, nn, nn );
        if all( sum( abs( xA - cyA * adiag ), 2 ) <= 2 * eps * sum( abs( xA ), 2 ) ),
            beta  = xB - alpha .* conj( yB );
            alpha = reshape( alpha, sz );
            if isreal( y ),
                cvx_optval = alpha .* square( y ) + reshape( beta, sz ) .* y;
            elseif all( abs( beta ) <= 2 * eps * abs( xB ) ),
                cvx_optval = alpha .* square_abs( y );
            else
                error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form(s): product is not real.\n' ) );
            end
        else
            error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form(s): not a square.\n' ) );
        end

    case 5,

        % posynomial .* posynomial
        cvx_optval = exp( log( xt ) + log( yt ) );

    otherwise,

        error( 'Shouldn''t be here.' );

    end

    %
    % Store the results
    %

    if nv == 1,
        z = cvx_optval;
    else
        z = cvx_subsasgn( z, t, cvx_optval );
    end

end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
