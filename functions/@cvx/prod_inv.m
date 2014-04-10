function y = prod_inv( varargin ) % x, dim, p

%PROD_INV   Internal cvx version.

if nargin < 3,
    varargin{3} = 1;
else
    p = varargin{3};
    if ~isnumeric(p) || ~isreal(p) || numel(p) ~= 1 || p <= 0 || isnan(p) || isinf(p),
        error( 'Third argument must be a positive scalar.' );
    end
end

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'concave' ; 'l_convex' ; 'l_concave' } );
    params.funcs = { @prod_inv_1, @prod_inv_2, @prod_inv_3, @prod_inv_3 };
    params.zero = 1;
    params.reduce = true;
    params.reverse = false;
    params.dimarg = 2;
    params.name = 'prod_inv';
end

try
    y = reduce_op( params, varargin{:} );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function y = prod_inv_1( x, p )
y = cvx( prod( cvx_constant( x ) ) .^ -p );

function y = prod_inv_2( x, p )
[nx,nv] = size( x ); %#ok
if nx == 1,
    y = cvx_recip( pos( x ) );
else
    y = [];
    cvx_begin
        epigraph variable y( 1, nv )
        geo_mean( [ x ; y ], 1, [ ones(nx,1) ; p ] ) >= 1; %#ok
        cvx_setnneg( y );
    cvx_end
end

function y = prod_inv_3( x, p )
y = exp( -p * sum( log( x ) ) );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
