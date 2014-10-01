function z = mtimes( x, y, oper )

%   Disciplined convex programming information for MTIMES:
%      True matrix multiplications Z = X * Y---that is, where neither X
%      nor Y is a scalar---require both multiplications and additions. 
%      For example, element (i,j) of Z is given by
%         X(i,1) * Y(1,j) + X(i,2) * Y(2,j) + ... + X(i,k) * Y(k,j)
%      Therefore, matrix multiplications must satisfy *both* the 
%      "no-product" rule of multiplication and the "same curvature" rule
%      of addition. See the help for CVX/TIMES and CVX/PLUS, 
%      respectively, for individual descriptions of these rules.
%   
%      An exception is made to these general rules for quadratic forms.
%      That is, two affine expressions may be multiplied together if the
%      result can immediately be verified as a convex quadratic form. 
%      For example, the construction
%         variable x(n)
%         x' * Q * x  <= 1;
%      would be permitted if Q is *constant* and positive semidefinite.
%   
%   Disciplined geometric programming information for TIMES:
%      As mentioned above, true matrix multiplies Z = X * Y require both
%      multiplications and additions. Since only log-convex terms can be
%      summed, both X and Y must be elementwise log-convex/affine.

if nargin < 3, oper = '*'; end

%
% Check sizes
%

sx = size( x );
sy = size( y );
if all( sx == 1 ) || all( sy == 1 ),
    z = times( x, y, oper );
    return
elseif length( sx ) > 2 || length( sy ) > 2,
    cvx_throw( 'Input arguments must be 2-D, or at least one input must be scalar.' );
elseif sx( 2 ) ~= sy( 1 ),
    cvx_throw( 'Inner matrix dimensions do not agree.' );
elseif sx( 2 ) == 1,
    z = bsxfun( @times, x, y, oper );
    return
end
sz = [ sx( 1 ), sy( 2 ) ];
nz = prod( sz );

%
% Check expression types
%

if cvx_isconstant( x ),
    x = sparse( cvx_constant( x ) );
    if cvx_isconstant( y ),
        % constant * constant
        y = sparse( cvx_constant( y ) );
        switch oper,
        case '\',  z = x \ y;
        case '/',  z = x / y;
        otherwise, z = x * y;
        end
        z = cvx( z );
    else
        % constant * everything
        z = y.basis_(:,reshape(1:prod(sy),sy)');
        z = reshape( z, [], sy(1) );
        switch oper,
        case '\', z = z / sparse( x.' );
        case '*', z = z * sparse( x.' );
        end
        z = reshape( z, [], nz );
        z = z(:,reshape(1:prod(sz),[sz(2),sz(1)])');
        z = cvx( sz, z );
    end
    if ~cvx_isvalid( z ),
        cvx_dcp_error( '*', 'mmult', x, y );
    end
    return
end
    
if cvx_isconstant( y ),
    % everything * constant
    y = sparse( cvx_constant( y ) );
    z = reshape( x.basis_, [], sx(2) );
    switch oper,
    case '/',  z = z / sparse( y );
    otherwise, z = z * sparse( y );
    end
    z = reshape( z, [], nz );
    z = cvx( sz, z );
    if ~cvx_isvalid( z ),
        cvx_dcp_error( '*', 'mmult', x, y );
    end
    return
end

vx = cvx_classify( x );
vy = cvx_classify( y );

persistent affnnc lvalid affine
if isempty( affnnc ),
    affnnc = cvx_remap( 'affine', 'n_concave', 'p_convex' );
    affine = cvx_remap( 'affine' );
    lvalid = cvx_remap( 'l_valid', 'nonnegative' );
end

if all(affnnc(vx)) && all(affnnc(vy)),
    if nz ~= 1,
        cvx_throw( 'Disciplined convex programming error:\n    Invalid quadratic form: not a scalar.\n' );
    end
    % quadratic form test 1: look for x' * y, x = D * conj( y ) + b
    % where D is diagonal psd or nsd, b is constant
    xA = x.basis_; 
    yA = y.basis_;
    mm = max( size( xA, 1 ), size( yA, 1 ) );
    xA( end + 1 : mm, : ) = 0;
    yA( end + 1 : mm, : ) = 0;
    xB  = xA( 2 : end, : );
    yB  = yA( 2 : end, : );
    cyB = conj( yB );
    alpha = sum( yB .* cyB, 1 );
    alpha = sum( bsxfun( @times, xB, yB ), 1 ) ./ max( alpha, realmin );
    if ~nnz( xB - bsxfun( @times, alpha, cyB ) > 2 * eps * sqrt( conj( xB ) .* xB ) ),
        if any( abs(imag(alpha)) > 2 * eps * abs(real(alpha)) ),
            cvx_throw( 'Disciplined convex programming error:\n    Invalid quadratic form: not real.\n' );
        end
        alpha = real( alpha );
        neg = any( alpha < 0 );
        if neg && ~all( alpha <= 0 ),
            cvx_throw( 'Disciplined convex programming error:\n    Invalid quadratic form: neither convex or concave.\n' );
        else
            offset = conj(x.basis_(1,:)-alpha.*y.basis_(1,:));
            if any( abs(imag(offset)) > 2*eps*abs(real(offset)) ),
                cvx_throw( 'Disciplined convex programming error:\n    Invalid quadratic form: not real.\n' );
            end
            z = y;
            if ~isreal( z ), z = abs(z); end
            z.basis_ = bsxfun( @times, sqrt(abs(alpha)), z.basis_ );
            z = sum_square( z );
            if neg, z.basis_ = -z.basis_; end
            offset = real(offset);
            if any(offset), z = z + offset * y; end
            return
        end
    end
    if all(affine(vx)) && all(affine(vy))
        dx = find( any( xA, 2 ) | any( yA, 2 ) );
        zb  = length( dx );
        xA  = xA( dx, : );
        yA  = yA( dx, : );
        P   = xA * yA.';
        Q   = xA * y.basis_(1,:).' + yA * x.basis_(1,:).';
        R   = x.basis_(1,:) * y.basis_(1,:).';
        if ~isreal( R ) || ~isreal( Q ) || ~isreal( P ),
            cvx_throw( 'Disciplined convex programming error:\n   Invalid quadratic form: not real.' );
        end
        xx = cvx( zb, sparse( dx, 1 : zb, 1 ) );
        [ z, success ] = quad_form( xx, P, Q, R );
        if ~success,
            cvx_throw( 'Disciplined convex programming error:\n    Invalid quadratic form: neither convex or concave.\n' );
        end
        return
    end
end

if all( lvalid(vx) ) && all( lvalid(vy) ),
    z = bsxfun( @times, x, reshape( y, [ 1, size(y) ] ) );
    z = sum( z, 2 );
    z = reshape( z, sz );
else
    cvx_dcp_error( '*', 'mmult', x, y );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
