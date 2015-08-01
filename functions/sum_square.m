function y = sum_square( X, dim, is_abs )

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

persistent P cname rname cmap rmap
if isempty( P ),
    cmap = cvx_remap( { ...
        { 'real', 'complex' } ; 'l_convex' ; ...
        { 'p_convex', 'n_concave', 'r_affine' } ; ...
        { 'p_convex', 'n_concave', 'affine' } } );
    cname = 'sum_square_abs';
    rmap = cvx_remap( { ...
        { 'real'  } ; 'l_convex' ; ...
        { 'p_convex', 'n_concave', 'r_affine' } } );
    rname = 'sum_square';
    P.funcs = { @ssqa_1, @ssqa_2, @ssqa_3, @ssqa_4 };
    P.zero = 0;
    P.constant = 1;
    P.reduce = true;
    P.reverse = false;
    P.dimarg = 2;
end
if nargin < 3 || ~is_abs,
    if nargin < 2, dim = []; end
    P.name = rname;
    P.map  = rmap;
else
    P.name = cname;
    P.map  = cmap;
end
y = cvx_reduce_op( P, X, dim );

function x = ssqa_1( x )
x = sum( x .* conj(x), 1 );

function x = ssqa_2( x )
x = sum( exp( 2 * log( x ), 1 ) );

function y = ssqa_3( x ) %#ok
[ nx, nv ] = size(x);
cvx_begin
    epigraph variable y( 1, nv ) nonnegative_
    { cvx_linearize(x), 0.5, y } == rotated_lorentz( [ nx, nv ], 1 ); %#ok
cvx_end

function y = ssqa_4( x )
y = ssqa_3( [ real(x) ; imag(x) ] );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
