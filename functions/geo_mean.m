function y = geo_mean( varargin )

%GEO_MEAN   Geometric mean.
%   Y=GEO_MEAN(X), where X is a vector, computes the geometrix mean of X. If any
%   of the elements of X are negative, then Y=-Inf. Otherwise, it is equivalent
%   to Y=PROD(X).^(1/LENGTH(X)). All elements must be real.
%
%   For matrices, GEO_MEAN(X) is a row vector containing the geometric means of
%   the columns. For N-D arrays, GEO_MEAN(X) is an array of the geometric means
%   taken along the first non-singleton dimension of X.
%
%   GEO_MEAN(X,DIM) takes the geometric mean along the dimension DIM of X.
%
%   GEO_MEAN(X,DIM,W), where W is a vector of nonnegative integers, computes a
%   weighted geometric mean Y = PROD(X.^W)^(1/SUM(W)). This is more efficient
%   than replicating the values of X W times. Note that W must be a vector,
%   even if X is a matrix, and its length must be the same as SIZE(X,DIM).
%
%   Disciplined convex programming information:
%       GEO_MEAN is concave  and nondecreasing; therefore, when used in CVX
%       specifications, its argument must be concave.

%
% Check arguments
%

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'concave' ; 'l_convex' ; 'l_concave' } );
    params.map = bsxfun( @and, params.map, ~cvx_remap( 'negative' ) );
    params.funcs = { @geo_mean_1, @geo_mean_2, @geo_mean_2, @geo_mean_2 };
    params.zero = 1;
    params.constant = 1;
    params.reduce = true;
    params.reverse = false;
    params.name = 'geo_mean';
    params.errargs = [];
    params.dimarg = [];
end

try
    [ sx, x, dim, w ] = cvx_get_dimension( 2, varargin );
    if ~isempty( w ),
        if ~( numel(w)==length(w) && isnumeric(w) && isreal(w) && all(w>=0) ),
            error( 'CVX:ArgError', 'Third argument must be a vector of nonnegative numbers.' );
        elseif ~any( w ),
            error( 'CVX:ArgError', 'The weight vector cannot be all zeros.')
        elseif numel( w ) ~= sx(dim),
            error( 'CVX:ArgError', 'Third argument must be a vector of length %d', sx(dim) );
        else
            w = w(:) / sum(w);
        end
    end
    y = cvx_reduce_op( params, x, dim, w );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function y = geo_mean_1( x, w )
[ nx, nv ] = size( x );
if isempty( w ), 
    w = 1 / nx;
else
    w = repmat( w, [ 1, nv ] ); 
end
y = prod( x .^ w, 1 );

function y = geo_mean_2( x, w ) %#ok
[ nx, nv ] = size( x );
cvx_begin
    hypograph variable y(1,nv);
    { cvx_accept_concave(x), y } == geo_mean_cone( [nx,nv], 1, w, 'func' ); %#ok
    cvx_setnneg(y);
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.% Copyright 2005-2014 CVX Research, Inc. 
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
