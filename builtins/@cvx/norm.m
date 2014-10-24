function cvx_optval = norm( x, p )

%   Disciplined convex programming information:
%       NORM is convex, except when P<1, so an error will result if
%       these non-convex "norms" are used within CVX expressions. NORM 
%       is nonmonotonic, so its input must be affine.

%
% Argument map
%

persistent remap1 remap2 remap3
if isempty( remap3 ),
    remap1 = cvx_remap( 'l_convex' );
    remap2 = cvx_remap( 'affine' );
    remap3 = cvx_remap( 'p_convex', 'n_concave', 'affine' );
end

%
% Check arguments
%

if nargin < 2 || isempty( p ),
    p = 2;
elseif ~isequal( p, 'fro' ) && ( ~isnumeric( p ) || ~isreal( p ) || p < 1 ),
    cvx_throw( 'Second argument must be a real number between 1 and Inf, or ''fro''.' );
elseif ndims( x ) > 2, %#ok
    cvx_throw( 'norm is not defined for N-D arrays.' );
end

[ m, n ] = size(x);
if m == 1 || n == 1 || isequal( p, 'fro' ),
    
    %
    % Vector norms
    %
    
    if isempty( x ),
        cvx_optval = cvx( 0 );
        return
    end
    if isequal( p, 'fro' ) || ~isreal(x) && p == 2
        p = 2;
        [ xR, x ] = bcompress( x );
        x = x .* sqrt(sum(xR.*conj(xR),2));
    else
        x = vec( x );
    end
    n = length( x );
    xc = cvx_classify( x );
    if ~all( remap3( xc ) ),
        cvx_throw( 'Disciplined convex programming error:\n    Cannot perform the operation norm( {%s}, %g )', cvx_class( x ), p );
    end
    if n == 1,
        cvx_optval = abs( x );
        return
    end
    switch p,
        case 1,
            cvx_optval = sum( abs( x ) );
        case Inf,
            cvx_optval = max( abs( x ) );
        otherwise,
            tt = remap1( xc );
            if all( tt ),
                cvx_optval = ( sum( x .^ p ) ) .^ ( 1 / p );
            else
                if nnz( tt ) > 1,
                    tt = tt ~= 0;
                    xx = cvx_subsref( x, tt );
                    xx = ( sum( xx .^ p ) ) .^ (1/p);
                    x  = [ cvx_subsref( x, ~tt ) ; xx ];
                end
                if ~all( remap2( xc ) ),
                    x = cvx_linearize( x );
                end
                if p == 2,
                    z = [];
                    cvx_begin
                        epigraph variable z nonnegative_
                        { x, z } == lorentz( n, [], ~isreal( x ) ); %#ok
                    cvx_end
                else
                    if isreal( x ),
                        cmode = 'abs';
                    else
                        cmode = 'cabs';
                    end
                    y = []; z = [];
                    cvx_begin
                        epigraph variable z nonnegative_
                        variable y( n )
                        { [ y, z*ones(n,1) ], x } == geo_mean_cone( [n,2], 2, [1/p,1-1/p], cmode ); %#ok
                        sum( y ) == z; %#ok
                    cvx_end
                end
            end
    end
    
else
    
    %
    % Matrix norms
    %
    
    xc = cvx_classify( x );
    if p == 2,
        map = remap2;
    else
        map = remap3;
    end
    if ~all( map( xc ) ),
        cvx_throw( 'Disciplined convex programming error:\n    Cannot perform the operation norm( {%s}, %g )\n   when the first argument is a matrix.', cvx_class( xt ), p );
    end
    switch p,
        case 1,
            cvx_optval = max( sum( abs( x ), 1 ), [], 2 );
        case Inf,
            cvx_optval = max( sum( abs( x ), 2 ), [], 1 );
        case 2,
            cvx_optval = sigma_max( x );
        otherwise,
            cvx_throw( 'The only matrix norms available are 1, 2, Inf, and ''fro''.' );
    end
    
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
