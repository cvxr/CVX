function y = geo_mean( x, dim, w )

%GEO_MEAN   Internal cvx version.

if nargin < 2 || isempty( dim ),
    dim = find( zx ~= 1, 1, 'first' );
    if isempty( dim ), dim = 1; end
elseif ~( isnumeric(dim) && numel(dim)==1 && isreal(dim) && dim > 0 && dim==floor(dim) ),
    error( 'Dimension argument must be a positive integer.' );
end

if nargin < 3,
    w = [];
elseif numel( w ) ~= length( w ) || ~isnumeric( w ) || ~isreal( w ) || any( w < 0 ) || any( w ~= floor( w ) ),
    error( 'Third argument must be a vector of nonnegative integers.' );
elseif length( w ) ~= size( x, dim ),
    error( 'Third argument must be a vector of length %d', nx );
end

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'concave' ; 'l_convex' ; 'l_concave' } );
    params.funcs = { @geo_mean_1, @geo_mean_2, @geo_mean_3, @geo_mean_3 };
    params.zero = 1;
    params.reduce = true;
    params.dimarg = 2;
    params.reverse = false;
    params.name = 'geo_mean';
end

try
    y = reduce_op( params, x, dim, w );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function y = geo_mean_1( x, w )
y = cvx( geo_mean( cvx_constant( x ), 1, w ) );

function y = geo_mean_2( x, w )
if nargin < 2, w = []; end
[nx,nv] = size(x);
y = [];
cvx_begin
    hypograph variable y(1,nv);
    { cvx_accept_concave(x), y } == geo_mean_cone( [nx,nv], 1,  w, 'func' ); %#ok
    cvx_setnneg(y);
cvx_end

function y = geo_mean_3( x, w )
nx = size(x);
if nx == 1,
    y = xt;
elseif nargin < 2 || isempty( w ),
    y = exp( sum( log( x ), 1 ) * ( 1 / nx ) );
else
    y = exp( ( w / sum( w ) ) * log( x ) );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
