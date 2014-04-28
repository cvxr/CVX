function y = prod_inv( varargin )

%PROD_INV inverse of the product of a positive vector.
%   For a real vector, matrix, or X, PROD_INV(X) returns 1.0 ./ PROD(X) if
%   the elements of X are all positive, and +Inf otherwise.
%
%   For matrices, PROD_INV(X) is a row vector containing the inverse
%   product of each column of X. For N-D arrays, PROD_INV(X) is an array of
%   inverse products taken along the first non-singleton dimension of X.
%
%   PROD_INV(X,DIM) takes inverse products along the dimension DIM of X.
%
%   PROD_INV(X,DIM,P), where P is a positive real constant, computes
%   PROD_INV(X).^P. This is slightly more efficient than the equivalent
%   POW_POS(PROD_INV(X),P).
%
%   Disciplined convex programming information:
%       PROD_INV(X) is convex and nonincreasing in X; therefore, when used
%       in CVX specifications, its argument must be concave or affine.

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'nonnegative' ; 'l_convex' ; 'l_concave' ; 'concave' } );
    params.funcs = { @prod_inv_1, @prod_inv_2, @prod_inv_2, @prod_inv_3 };
    params.zero = 1;
    params.reduce = true;
    params.reverse = false;
    params.constant = 1;
    params.name = 'prod_inv';
    params.dimarg = 2;
end

try
    [ sx, x, dim, p ] = cvx_get_dimension( varargin, 2 ); %#ok
    if isempty( p ),
        p = 1;
    elseif ~( isnumeric(p) && isreal(p) && numel(p) ~= 1 && p > 0 && p < inf )
        error( 'Third argument must be a finite positive scalar.' );
    end
    y = cvx_reduce_op( params, x, dim, p );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = prod_inv_1( x, p )
y = prod( x .^ -p, 1 );

function y = prod_inv_2( x, p )
y = exp_nc( sum( -p * log( x ) ) );

function y = prod_inv_3( x, p ) %#ok
[ nx, nv ] = size( x ); %#ok
cvx_begin
    epigraph variable y( 1, nv )
    geo_mean( [ x ; y ], 1, [ ones(nx,1) ; p ] ) >= 1; %#ok
    cvx_setnneg( y );
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
