function z = mtimes( x, y )
error( cvx_verify( x, y ) );

%
% Quick exit for constant case
%

xc = cvx_isconstant( x );
yc = cvx_isconstant( y );
if xc & yc,
    z = builtin( 'mtimes', cvx_constant( x ), cvx_constant( y ) );
    return
end

%
% Check sizes
%

sx = size( x );
sy = size( y );
if all( sx == 1 ) | all( sy == 1 ),
    z = pluslikeoper( x, y, 'times' );
    return
elseif length( sx ) > 2 | length( sy ) > 2,
    error( 'Input arguments must be 2-D.' );
elseif sx( 2 ) ~= sy( 1 ),
    error( 'Inner matrix dimensions must agree.' );
else,
    sz = [ sx( 1 ), sy( 2 ) ];
end

[ prob, x, y ] = cvx_operate( [], x, y );
    
%
% Quadratic forms
%

if ~xc & ~yc,
    
    %
    % Both terms must be affine, and the result must be a scalar
    %
    
    if ~cvx_isaffine( x ) | ~cvx_isaffine( y ),
        error( sprintf( 'Disciplined convex programming error:\n    Multiplication of non-affine expressions is forbidden.\n' ) );
    elseif any( sz ~= 1 ),
        error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form: must be a scalar.' ) );
    end
    
    %
    % Separate the quadratic and non-quadratic terms
    %

    qq = cvx_isconstant( x, true )' | cvx_isconstant( y, true );
    if any( qq ),
        z = 0;
        qt.type = '()';
        for k = 1 : 2,
            qt.subs = { qq };
            xx = subref( x, qt.subs, 0 );
            yy = subref( y, qt.subs, 0 );
            z  = z + mtimes( xx( : ).' * yy( : ) );
            qq = ~qq;
        end
        return
    end
    
    %
    % Quadratic form test 1: See if x == a conj( y ) + b for some real a, b,
    % so that the quadratic form involves a simple squaring (or sum of squares)
    %

    mm    = dimension( prob ) - 1;
    xA    = x.basis_( :, 2 : end );
    yA    = y.basis_( :, 2 : end );
    if size( xA, 2 ) < mm, xA( end, mm ) = 0; end
    if size( yA, 2 ) < mm, yA( end, mm ) = 0; end
    cyA   = conj( yA );
    alpha = sum( sum( real( xA .* yA ) ) ) ./ max( sum( sum( cyA .* yA ) ), realmin );
    if sum( sum( abs( xA - alpha * cyA ) ) ) <= 2 * eps * sum( sum( abs( xA ) ) ),
        beta = x.basis_( :, 1 ) - alpha * conj( y.basis_( :, 1 ) );
        if isreal( y ) & isreal( beta ),
            z = alpha * sum_square( y ) + beta .' * y( : );
        elseif all( abs( beta ) <= 2 * eps * abs( x.basis_( :, 1 ) ) ),
            z = alpha * sum_square_abs( y );
        else,
            error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form: product is not real.\n' ) );
        end
        return
    end

    %
    % Quadratic form test 2: Extract the quadratic coefficient matrix
    % and test it for semidefiniteness
    %
    
    dx = find( any( xA, 1 ) | any( yA, 1 ) );
    zb = length( dx );
    xA = xA( :, dx )';
    yA = yA( :, dx )';
    xB = x.basis_( :, 1 );
    yB = y.basis_( :, 1 );
    P  = xA * yA.';
    Q  = xA * yB + yA * xB;
    R  = xB.' * yB;
    P  = 0.5 * ( P + P.' );
    if ~isreal( R ) | ~isreal( Q ) | ~isreal( P ),
        error( sprintf( 'Disciplined convex programming error:\n   Invalid quadratic form: product is complex.' ) );
    else,
        xx = cvx( prob, zb, sparse( 1 : zb, dx + 1, 1 ) );
        [ z, success ] = quad_form( xx, P );
        if ~success,
            error( sprintf( 'Disciplined convex programming error:\n   Invalid quadratic form: neither convex nor concave.' ) );
        end
        z = z + Q' * xx + R;
    end
    
else,
    
    if xc,
        bx = cvx_constant( x );
        by = cvx_basis( y );
    else,
        bx = cvx_constant( y ).';
        by = cvx_basis( x.' );
        sy = [ sx( 2 ), sx( 1 ) ];
        sz = [ sz( 2 ), sz( 1 ) ];
    end
    nv = size( by, 2 );
    by = reshape( by, [ sy( 1 ), sy( 2 ) * nv ] );
    z  = reshape( bx * by, [ prod( sz ), nv ] );
    z  = cvx( prob, sz, z );
    if ~xc,
        z = z.';
    end
    v = cvx_vexity( z );
    if any( isnan( v( : ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n   Illegal affine combination of convex and/or concave terms detected.' ) );
    end
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
