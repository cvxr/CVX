function cvx_optval = norm( x, p )
error( nargchk( 1, 2, nargin ) );

%NORM   Matrix or vector norm.
%
%   CVX supports all of the  matrix and vector norms that MATLAB itself
%   supports, with two caveats. First, norm(X,-Inf) does not compute a 
%   true norm, so CVX does not support this mode. Second, for norm(X,p),
%   1<p<+Inf, CVX approximates p with a rational value as computed by
%   the RAT() function. The larger the numerator of this ratio, the more
%   complex the computations that CVX initiates. See the CVX user guide
%   under "rational powers" for more details.
%
%   For vectors...
%     NORM(V,P)           = sum(abs(V).^P)^(1/P)
%     NORM(V)             = norm(V,2).
%     NORM(V,inf)         = max(abs(V)).
%     NORM(V,-inf)        = min(abs(V)).
%   Note that NORM(V,P) for -Inf<=P<1 is NOT, in fact, a valid norm. MATLAB
%   allows them nonetheless, so CVX does as well---but for numeric vectors
%   ONLY. Such "norms" cannot be used within disciplined convex programs.
%
%   For matrices...
%     NORM(X)             = max(svd(X)).
%     NORM(X,2)           = norm(X,2).
%     NORM(X,1)           = max(sum(abs(X))).
%     NORM(X,Inf)         = max(sum(abs(X'))).
%     NORM(X,'fro')       = norm(X(:),2).
%   NORM(X,P) is not implemented for matrices for any other values of P.
%
%   Disciplined convex programming information:
%       NORM is convex, except when P<1, so an error will result if these
%       non-convex "norms" are used within CVX expressions. NORM is
%       nonmonotonic, so its input must be affine.

%
% Check arguments
%

if nargin < 2,
    p = 2;
elseif ~isequal( p, 'fro' ) & ( ~isnumeric( p ) | ~isreal( p ) | p < 1 ),
    error( 'Second argument must be a real number between 1 and Inf, or ''fro''.' );
end
if ndims( x ) > 2,
    error( 'norm is not defined for N-D arrays.' );
end

[ m, n ] = size( x );
if m == 0 | n == 0,
    
    %
    % Empty matrices
    %
    
    cvx_optval = 0;
    
elseif m == 1 & n == 1,
    
    %
    % Scalars
    %
    
    cvx_optval = abs( x );
    
elseif m == 1 | n == 1 | isequal( p, 'fro' ),
    
    %
    % Vector norms
    %

    if isequal( p, 'fro' ), p = 2; end
    x = svec( x, p );
    persistent remap1 remap2 remap3
    if isempty( remap2 ),
        remap1 = cvx_remap( 'log-convex' );
        remap2 = cvx_remap( 'affine', 'log-convex' );
        remap3 = cvx_remap( 'log-affine' );
    end
    xc = cvx_classify( x );
    if ~all( remap2( xc ) ) & ~( isequal( p, 1 )& ~all( remap3( xc ) ) ),
        error( sprintf( 'Disciplined convex programming error:\n    Cannot perform the operation norm( {%s}, %g )', cvx_class( xt ), p ) );
    end
    switch p,
        case 1,
            cvx_optval = sum( abs( x ) );
        case Inf,
            cvx_optval = max( abs( x ) );
        otherwise,
            tt = remap1( xc );
            if all( tt ),
                cvx_optval = ( sum( x .^ p ) ) .^ (1/p);
            else
                if nnz( tt ) > 1,
                    tt = tt ~= 0;
                    xx = cvx_subsref( x, tt );
                    xx = ( sum( xx .^ p ) ) .^ (1/p);
                    x  = [ cvx_subsref( x, ~tt ) ; cvx_accept_convex( xx ) ];
                end
                n = length( x );
                if p == 2,
                    cvx_begin
                        epigraph variable z
                        { x, z } == lorentz( n, [], ~isreal( x ) );
                    cvx_end
                else
                    map = cvx_geomean_map( p, true );
                    cvx_begin
                        epigraph variable z
                        if rem( p, 2 ) == 0,
                            p = p * 0.5;
                            x_abs = quad_over_lin( x, z, 0 );
                        else
                            x_abs = abs( x );
                        end
                        variable y( n )
                        geomean( [ y, z*ones(n,1) ], 2, map, true ) >= x_abs;
                        sum( y ) == z;
                    cvx_end
                end
            end
    end
    
else
    
    %
    % Matrix norms
    %
    
    if ~cvx_isaffine( x ),
        error( sprintf( 'Disciplined convex programming error:\n    Cannot perform the operation norm( {%s}, %g )\n   when the first argument is a matrix.', cvx_class( xt ), p ) );
    end
    switch p,
        case 1,
            cvx_optval = max( sum( abs( x ), 1 ), [], 2 );
        case Inf,
            cvx_optval = max( sum( abs( x ), 2 ), [], 1 );
        case 2,
            cvx_optval = sigma_max( x );
        otherwise,
            error( 'The only matrix norms available are 1, 2, Inf, and ''fro''.' );
    end
    
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
