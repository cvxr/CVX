function y = sum_square( varargin )

%SUM_SQUARE   Sum of squares.
%   For vectors, SUM_SQUARE(X) is the sum of the squares of the elements of
%   the vector; i.e., SUM(X.^2).
%
%   For matrices, SUM_SQUARE(X) is a row vector containing the application
%   of SUM_SQUARE to each column. For N-D arrays, the SUM_SQUARE operation
%   is applied to the first non-singleton dimension of X.
%
%   SUM_SQUARE(X,DIM) takes the sum along the dimension DIM of X.
%
%   Disciplined convex programming information:
%       If X is real, then SUM_SQUARE(X,...) is convex and nonmonotonic in
%       X. If X is complex, then SUM_SQUARE(X,...) is neither convex nor
%       concave. Thus, when used in CVX expressions, X must be affine. DIM
%       must be constant.

persistent params
if isempty( params ),
    params.map = cvx_remap( { 'real' ; 'l_convex' ; ...
         { 'p_convex', 'n_concave', 'r_affine' } } );
    params.funcs = { @ssq_1, @ssq_2, @ssq_3 };
    params.zero = 0;
    params.constant = 1;
    params.reduce = true;
    params.reverse = false;
    params.name = 'sum_square';
    params.dimarg = 2;
end

try
    y = cvx_reduce_op( params, varargin{:} );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end

function x = ssq_1( x )
x = sum( x .^ 2, 1 );

function x = ssq_2( x )
x = sum( exp_nc( 2 * log( x ), 1 ) );

function y = ssq_3( x ) %#ok
[ nx, nv ] = size(x);
cvx_begin
    epigraph variable y( 1, nv ) nonnegative_
    { linearize_abs(x), y, 1 } == rotated_lorentz( [ nx, nv ], 1 ); %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
