function z = pluslikeoper( x, y, oper, cheat )
if nargin < 4, cheat = false; end
error( cvx_verify( x, y ) );

%
% Quick exit for constant case
%

xc = cvx_isconstant( x );
yc = cvx_isconstant( y );
if xc & yc,
    z = builtin( oper, cvx_constant( x ), cvx_constant( y ) );
    return
end

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
else,
    sz = sx;
end

%
% Quadratic form test
%

if ~xc & ~yc & isequal( oper, 'times' ),
    
    %
    % Separate into quadratic and non-quadratic cases
    %
    
    xN = cvx_isconstant( x, true );
    yN = cvx_isconstant( y, true );
    qq = xN | yN;
    if any( qq ),
        z = zeros( sz );
        qt.type = '()';
        for k = 1 : 2,
            qt.subs = { qq };
            if xs, xx = x; else, xx = subsref( x, qt ); end
            if ys, yy = y; else, yy = subsref( y, qt ); end
            z = subsasgn( z, qt, pluslikeoper( xx, yy, oper, cheat ) );
            qq = ~qq;
        end
        return
    end
    
    %
    % Terms must be affine
    %
    
    if ~cvx_isaffine( x ) | ~cvx_isaffine( y ),
        error( sprintf( 'Disciplined convex programming error:\n    Multiplication of non-affine expressions is forbidden.\n' ) );
    end
            
    %
    % Quadratic form test: See if x == a conj( y ) + b for some real a, b, so that
    % the quadratic form involves a simple squaring
    %
  
    [ prob, x, y ] = cvx_operate( [], x, y );
    nn    = prod( sz );
    mm    = dimension( prob ) - 1;
    xA    = x.basis_( :, 2 : end );
    yA    = y.basis_( :, 2 : end );
    if size( xA, 2 ) < mm, xA( end, mm ) = 0; end
    if size( yA, 2 ) < mm, yA( end, mm ) = 0; end
    cyA   = conj( yA );
    alpha = sum( real( xA .* yA ), 2 ) ./ max( sum( cyA .* yA, 2 ), realmin );
    adiag = sparse( 1 : nn, 1 : nn, alpha, nn, nn );
    if all( sum( abs( xA - adiag * cyA ), 1 ) <= 2 * eps * sum( abs( xA ), 1 ) ),
        beta  = x.basis_( :, 1 ) - alpha .* conj( y.basis_( :, 1 ) );
        alpha = reshape( alpha, sz );
        if isreal( y ),
            z = alpha .* square( y ) + reshape( beta, sz ) .* y;
        elseif all( abs( beta ) <= 2 * eps * abs( x.basis_( 1, : ) ) ),
            z = alpha .* square_abs( y );
        else,
            error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form(s): product is not real.\n' ) );
        end
    else,
        error( sprintf( 'Disciplined convex programming error:\n    Invalid quadratic form(s): not a square.\n' ) );
    end
    
    return
    
end

%
% Check vexity
%

if ~cheat,
    switch oper,
        case 'plus',
            if ~xc & ~yc,
                temp = cvx_vexity( x ) .* cvx_vexity( y );
                if any( temp( : ) < 0 ),
                    error( sprintf( 'Disciplined convex programming error:\n   Addition of convex and concave terms is forbidden.' ) );
                end
            end
        case 'minus',
            if ~xc & ~yc,
                temp = cvx_vexity( x ) .* cvx_vexity( y );
                if any( temp( : ) > 0 ),
                    error( sprintf( 'Disciplined convex programming error:\n   Differences between convex terms or concave terms are forbidden.' ) );
                end
            end
        case 'ldivide',
            if ~xc,
                error( sprintf( 'Disciplined convex programming error:\n   Division by a non-constant expression is forbidden.' ) );
            elseif any( x.basis_( 1, : ) == 0 ),
                error( 'Division by zero.' );
            end
        case 'rdivide',
            if ~yc,
                error( sprintf( 'Disciplined convex programming error:\n   Division by a non-constant expression is forbidden.' ) );
            elseif any( y.basis_( 1, : ) == 0 ),
                error( 'Division by zero.' );
            end
    end
end

%
% Apply operation, stretching basis matrices as needed
%

[ prob, x, y ] = cvx_operate( [], x, y );
bx = cvx_basis( x );
by = cvx_basis( y );
switch oper,
    case { 'times', 'rdivide', 'ldivide' },
        nn = prod( sz );
        if xc,
            bx = bx( :, 1 );
            if ~isequal( oper, 'times' ), bx = 1.0 ./ bx; end
            if ~xs, bx = sparse( 1 : nn, 1 : nn, bx, nn, nn ); end
            bz = bx * by;
        else, % yc
            by = by( :, 1 );
            if ~isequal( oper, 'times' ), by = 1.0 ./ by; end
            if ~xs, by = sparse( 1 : nn, 1 : nn, by, nn, nn ); end
            bz = by * bx;
        end
    otherwise,
        dx = size( bx );
        dy = size( by );
        if dx( 2 ) < dy( 2 ),
            bx( end, dy( 2 ) ) = 0;
        elseif dx( 2 ) > dy( 2 ),
            by( end, dx( 2 ) ) = 0;
        end
        if xs & ~ys,
            bx = bx( ones( 1, dy( 1 ) ), : );
        elseif ys & ~xs,
            by = by( ones( 1, dx( 1 ) ), : );
        end
        if isequal( oper, 'minus' ), by = -by; end
        bz = bx + by;
end

%
% Construct result
%

z = cvx( prob, sz, bz );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
