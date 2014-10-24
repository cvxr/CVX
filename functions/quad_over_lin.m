function z = quad_over_lin( varargin )

%QUAD_OVER_LIN Sum of squares over linear.
%   Z=QUAD_OVER_LIN(X,Y), where X is a vector and Y is a scalar, is equal to
%   SUM(ABS(X).^2)./Y if Y is positive, and +Inf otherwise. Y must be real.
%
%   If X is a matrix, QUAD_OVER_LIN(X,Y) is a row vector containing the values
%   of QUAD_OVER_LIN applied to each column. If X is an N-D array, the operation
%   is applied to the first non-singleton dimension of X.
%
%   QUAD_OVER_LIN(X,Y,DIM) takes the sum along the dimension DIM of X.
%   A special value of DIM == 0 is accepted here, which is automatically
%   replaced with DIM == NDIMS(X) + 1. This has the effect of eliminating
%   the sum; thus QUAD_OVER_LIN( X, Y, NDIMS(X) + 1 ) = ABS( X ).^2 ./ Y.
%
%   In all cases, Y must be compatible in the same sense as ./ with the squared
%   sum; that is, Y must be a scalar or the same size as SUM(ABS(X).^2,DIM).
%
%   Disciplined convex programming information:
%       QUAD_OVER_LIN is convex, nonmontonic in X, and nonincreasing in Y.
%       Thus when used with CVX expressions, X must be convex (or affine)
%       and Y must be concave (or affine).

%
% Check arguments
%

persistent P
if isempty( P ),
    P.map = cvx_remap( ...
        { { 'any' }, { 'nonpositive' } }, ...
        { { 'constant' }, { 'positive' } }, ...
        { { 'l_convex' },  { 'l_concave' } }, ...
        { { 'l_concave' }, { 'l_convex' } }, ...
        { { 'r_affine', 'p_convex', 'n_concave' }, { 'positive' } }, ...
        { { 'r_affine', 'p_convex', 'n_concave' }, { 'concave' } }, ...
        { { 'affine', 'p_convex', 'n_concave' }, { 'positive' } }, ...
        { { 'affine', 'p_convex', 'n_concave' }, { 'concave' } }, ...
        [ 0, 1, 2, 2, 3, 4, 5, 6 ] );
    P.funcs = { @qol_cnst, @qol_log, @qol_sqr, @qol_lin, @qol_sqa, @qol_lin, @qol_cpx };
    P.constant = 1;
    P.name = 'quad_over_lin';
end
[ sx, x, y, dim ] = cvx_get_dimension( varargin, 3, 'zero', true );
if sx(dim) > 1, x = norms( x, 2, dim ); end
z = cvx_binary_op( P, x, y );

function z = qol_cnst( x, y )
z = x .^ 2 / y;

function z = qol_log( x, y )
z = exp( 2 * log( x ) - log( y ) );

function z = qol_sqr( x, y )
z = square( x ) ./ y;

function z = qol_lin( x, y ) %#ok
sz = max( size(x), size(y) );
cvx_begin
    epigraph variable z( sz ) nonnegative_
    { cvx_linearize(x), cvx_linearize(y), 0.5 * z } == rotated_lorentz( sz, 0 ); %#ok
cvx_end

function z = qol_sqa( x, y )
z = square( abs( x ) ) ./ y;

function z = qol_cpx( x, y ) %#ok
sz = max( size(x), size(y) );
cvx_begin
    epigraph variable z( sz ) nonnegative_
    { cvx_linearize(x), cvx_linearize(y), 0.5 * z } == rotated_complex_lorentz( sz, 0 ); %#ok
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
