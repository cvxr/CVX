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
    
    qx = cvx_isconstant( x, true )';
    qy = cvx_isconstant( y, true );
    if any( qx | qy ),
        z = 0;
        qt.type = '()';
        for k = 0 : 3,
            qt.subs = { qx & qy };
            if any( qt.subs{1} ),
                z = z + subsref( x, qt ) * subsref( y, qt );
            end
            if rem( k, 2 ) == 1, 
                qy = ~qy; 
            end
            qx = ~qx;
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
    
    %
    % Eliminate any zero rows from the left-hand matrix
    %
    
    if xc,
        bx = cvx_constant( x );
        by = cvx_basis( y );
    else,
        bx = cvx_constant( y ).';
        by = cvx_basis( x );
        sy = sx;
    end
    tx = any( bx, 2 );
    if nnz( tx ) ~= numel( tx ),
        bx = bx( tx, : );
        tx = full( find( tx ) );
    else,
        tx = [];
    end
    
    %
    % Determine the 3-D indices for the right-hand basis,
    % and reshape to sy( 1 ) x ( sy( 2 ) * nv )
    %

    nv = size( by, 2 );
    [ yi, yb, by ] = find( by );
    yj = floor( ( yi - 1 ) / sy( 1 ) );
    yi = yi - sy( 1 ) * yj;
    yj = yj + 1;
    if xc,
        yj = yj + ( yb - 1 ) * sy( 2 );
    else,
        yb = yi + ( yb - 1 ) * sy( 1 );
        yi = yj; yj = yb;
        sy = [ sy( 2 ), sy( 1 ) ];
    end
    clear yij yb
    
    %
    % Eliminate any zero columns from this reshaped right-hand matrix
    %
    
    ty  = sparse( yj, 1, 1, nv * sy( 2 ), 1 );
    nv2 = nnz( ty );
    if nv2 ~= numel( ty ),
        ty  = find( ty );
        sv  = sparse( ty, 1, 1 : nv2 );
        yj  = full( sv( yj ) );
        clear sv
    else,
        ty = [];
    end
    
    %
    % Perform the multiplication and re-expand
    %
    
    by = sparse( yi, yj, by, sy( 1 ), nv2 );
    clear yi yj
    z = bx * by;
    clear bx by
    [ zi, zj, z ] = find( z );
    if ~isempty( tx ), zi = tx( zi ); end; clear tx
    if ~isempty( ty ), zj = ty( zj ); end; clear ty
    
    %
    % Reshape back to ( sz( 1 ) * sz( 2 ) ) x nv
    %
    
    zb = floor( ( zj - 1 ) / sy( 2 ) );
    zj = zj - zb * sy( 2 );
    zb = zb + 1;
    if xc,
        zi = zi( : ) + sz( 1 ) * ( zj( : ) - 1 );
    else,
        zi = zj( : ) + sz( 1 ) * ( zi( : ) - 1 );
    end
    
    %
    % Construct the result
    %
    
    z = sparse( zi, zb, z, prod( sz ), nv );
    z = cvx( prob, sz, z );
    
    %
    % Check that the sums are legal
    %
    
    v = cvx_vexity( z );
    if any( isnan( v( : ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n   Illegal affine combination of convex and/or concave terms detected.' ) );
    end
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
