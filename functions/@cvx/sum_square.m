function y = sum_square( varargin ) % x, dim

%SUM_SQUARE   Internal cvx version.

persistent params
if isempty( params ),
    params.map = cvx_remap( { ...
        'real' ; ...
        'l_convex' ; ...
        { 'p_convex', 'n_concave', 'affine' } } );
    params.funcs = { @ssq_1, @ssq_2, @ssq_3 };
    params.zero = 0;
    params.reduce = true;
    params.reverse = false;
    params.dimarg = 2;
    params.name = 'sum_square';
end

try
    y = reduce_op( params, varargin{:} );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function x = ssq_1( x )
y = sum( x .^ 2 );

function x = ssq_2( x )
x = exp( 2 * log( x ) );
if x.size_ > 1, x = sum( x ); end

function y = ssq_3( x )
[nx,nv] = size(x);
y = [];
x = cvx_accept_cvxccv( x );
cvx_begin
    epigraph variable y( 1, nv )
    { x, y, 1 } == rotated_lorentz( [ nx, nv ], 1, ~isreal( x ) ); %#ok
    cvx_setnneg(y);
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
