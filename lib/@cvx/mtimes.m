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
    
    if xc,
        bx = cvx_constant( x );
        if sy(2) ~= 1, bx = kron( speye(sy(2)), bx ); end
        by = cvx_basis( y );
    else,
        bx = cvx_constant(y).';
        if sx(1) ~= 1, bx = kron( bx, speye(sx(1)) ); end
        by = cvx_basis( x );
    end
    nc = 8 * ( 1 + ~( isreal( bx ) & isreal( by ) ) );
    nc = max( floor( 128 * 1024 * 1024 / nc / size(by,1) ), 1 );
    if nc >= size(by,2),
        bz = bx * by;
    else,
        bz = {};
        tt = any( by, 1 );
        if nnz( tt ) > 0.5 * nnz( tt ), 
            tt( : ) = true; 
        end
        tv = find( tt );
        nt = length( tv );
        for k = rem( nt, nc ) : nc : nt,
            ttv = tv(max(k-nc,0)+1:k);
            bz{end+1} = bx * by( :, ttv );
        end
        bz = horzcat( bz{:} );
        if ~all( tt ),
            if issparse( bz ),
                [ zi, zj, zv ] = find( bz );
                bz = sparse( zi, tv(zj), zv, size(bz,1), size(by,2) ); 
            else,
                bz = cvx_subasgn( zeros( size(bz,1), size(by,2) ), ':', tt, bz );
            end
        end
    end
    z = cvx( prob, sz, bz );
    
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
