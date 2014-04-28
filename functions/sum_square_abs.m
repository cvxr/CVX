function y = sum_square_abs( varargin )

%SUM_SQUARE   Sum of squares.
%   For vectors, SUM_SQUARE_ABS(X) is the sum of the squares of the
%   absolute values of the elements of the vector; i.e., SUM(ABS(X).^2).
%
%   For matrices, SUM_SQUARE_ABS(X) is a row vector containing the
%   application of SUM_SQUARE_ABS to each column. For N-D arrays, the 
%   SUM_SQUARE_ABS operation is applied to the first non-singleton 
%   dimension of X.
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
    params.map = cvx_remap( { ...
        { 'real', 'complex' } ; 'l_convex' ; 'c_affine' ;
        { 'p_convex', 'n_concave', 'r_affine' } } );
    params.funcs = { @ssqa_1, @ssqa_2, @ssqa_3, @ssqa_3 };
    params.zero = 0;
    params.constant = 1;
    params.reduce = true;
    params.reverse = false;
    params.name = 'sum_square_abs';
    params.dimarg = 2;
end

try
    y = cvx_reduce_op( params, varargin{:} );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); a
    else rethrow( exc ); end
end

function x = ssqa_1( x )
x = sum( x .* conj(x), 1 );

function x = ssqa_2( x )
x = sum( exp_nc( 2 * log( x ), 1 ) );

function y = ssqa_3( x ) %#ok
[ nx, nv ] = size(x);
cvx_begin
    epigraph variable y( 1, nv ) nonnegative_
    { linearize_abs(x), y, 1 } == rotated_lorentz( [ nx, nv ], 1, ~isreal(x) ); %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.