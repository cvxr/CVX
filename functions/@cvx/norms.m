function y = norms( x, p, dim )

%NORMS   Internal cvx version.

if nargin < 2 || isempty(p),
    p = 2;
elseif ~isnumeric( p ) || numel( p ) ~= 1 || ~isreal( p ) || isnan( p ) || p < 1,
    error( 'Second argument must be a scalar between 1 and +Inf, inclusive.' );
end

try
    if size( x, dim ) == 0,
        y = abs( x );
        return
    else
        switch p,
            case 0,
                y = abs( x );
                return
            case 1,
                y = sum( abs( x ), dim );
                return
            case Inf,
                y = max( abs( x ), [], dim );
                return
        end
    end
catch
    cvx_dcp_error( x, p, 'norms' );
end

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'constant' ; 'l_convex' ; { 'p_convex', 'p_concave', 'affine' } } );
    params.funcs = { @norms_1, @norms_2, @norms_3 };
    params.zero = 0;
    params.reduce = true;
    params.reverse = false;
    params.dimarg = 3;
    params.fname = 'norms';
end

try
    y = reduce_op( params, x, p, dim );
catch exc
    if isequal( exc.identifier, 'CVX:DCPError' ), throw( exc ); 
    else rethrow( exc ); end
end

function y = norms_1( x, p )
y = cvx( sum( abs( cvx_constant( x ) ) .^ p ) .^ (1/p) );

function y = norms_2( x, p )
y = sum( abs( x ) .^ p ) .^ (1/p);

function y = norms_3( x, p )
[nx,nv] = size(x);
y = [];
x = cvx_accept_cvxccv( x );
if p == 2,
    cvx_begin
        epigraph variable y( 1, nv )
        { x, y } == lorentz( [ nx, nv ], 1, ~isreal( x ) ); %#ok
        cvx_setnneg(y);
    cvx_end
else
    z = []; y = [];
    cvx_begin
        variable z( nx, nv )
        epigraph variable y( 1, nv )
        if isreal(x), cmode = 'abs'; else cmode = 'cabs'; end
        { cat( 3, z, repmat(y,[nx,1]) ), x } ...
            == geo_mean_cone( [nx,nv,2], 3, [1/p,1-1/p], cmode ); %#ok
        sum( z ) == y; %#ok
        cvx_setnneg(y);
    cvx_end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
