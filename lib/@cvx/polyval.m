function y = polyval( p, x )

%POLYVAL    Evaluate polynomial.
%
%   Y = POLYVAL(P,X) returns the value of a polynomial P evaluated at X. P
%   is a vector of length N+1 whose elements are the coefficients of the
%   polynomial in descending powers.
% 
%       Y = P(1)*X^N + P(2)*X^(N-1) + ... + P(N)*X + P(N+1)
% 
%   If X is a matrix or vector, the polynomial is evaluated at all
%   points in X.
%
%   POLYVAL can be used with CVX variables in two ways:
%   --- If P is constant and X is a variable, then POLYVAL(P,X) represents
%       a polynomial function of X. The polynomial described by P must be
%       affine, convex, or concave. If LENGTH(P) > 2, then X must be real
%       and affine; otherwise, X may have any curvature.
%   --- If X is a constant and P is a variable, then POLYVAL(P,X) performs
%       a simple linear combination of the elements of P.
%   POLYVAL cannot accept variables for P and X simultaneously.
%
%   Disciplined convex programming information:
%       Refer to the description above for more information.

sp = size( p );
if isempty( p ),
    p = zeros( 1, 0 );
elseif length( sp ) > 2 | ~any( sp == 1 ),
    error( 'First argument must be a vector.' );
end
n = length( p );
sx = size(x);

if cvx_isconstant( p ),
    
    if ~isreal( p ),
        error( 'Polynomial must be real.' );
    end
    p = cvx_constant( p );
    if cvx_isconstant( x ),
        y = polyval( p, cvx_constant( x ) );
        return
    end
    if any( isinf( p ) | isnan( p ) ),
        error( 'Inf and NaN not accepted here.' );
    end
    ndxs = find( p );
    if isempty( ndxs ),
        y = zeros( sx );
        return
    else
        for k = ndxs(:)',
            pt = p( k : end );
            if ~any( isinf( pt ./ pt(1) ) ),
                p = pt;
                break;
            end
        end
    end
    n = length( p );
    switch min( n, 4 + rem( n, 2 ) ),
        case 1,
            % Constant
            y = p(1) * ones(sx);
        case 2,
            % Affine
            y = p(1) * x + p(2);
        case 3,
            % Quadratic
            b2a = p(2) ./ ( 2 * p(1) );
            y = p(1) * square( x + b2a ) + ( p(3) - b2a );
        case 4,
            % Odd powers
            error( 'Polynomial must be affine, convex, or concave.' );
        case 5,
            % Even powers
            pd = roots( [n-2:-1:1].*[n-1:-1:2].*reshape(p(1:end-2),1,n-2) );
            pd = pd( imag(pd) == 0 );
            if ~isempty( pd ),
                pr = diff( [ 0, find(diff(sort(pd))~=0), length(pd) ] );
                if any( rem( pr, 2 ) ),
                    error( 'Polynomial must be affine, convex, or concave.' );
                end
            end
            y = polyenv( p, x );
    end
    
elseif cvx_isconstant( x ),

    [ ii, jj, vv ] = find( p );
    jj = ii + jj - 1;
    nv = length(vv);
    nx = prod( sx );
    y = reshape( x, nx, 1 ) * ones( 1, nv );
    y = y .^ ( ones( nx, 1 ) * (n-jj(:))' );
    y = reshape( y * reshape( vv, nv, 1 ), sx );
    
else
    
    error( 'At least one of the arguments must be constant.' );
    
end
