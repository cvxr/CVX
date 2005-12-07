function z = quad_over_lin( x, y, dim )
error( nargchk( 2, 3, nargin ) );

%QUAD_OVER_LIN   quadratic over linear.
%
%   For real vectors X, QUAD_OVER_LIN(X,Y) is the sum of the squares of the
%   elements of X, divided by Y; i.e., SUM( X .^ 2 ) / Y.
%
%   For complex vectors X, QUAD_OVER_LIN(X,Y) is the sum of the squares of
%   the magnitudes of X, divided by Y; i.e., SUM( ABS( X ) .^ 2 ) / Y.
%
%   For matrices, QUAD_OVER_LIN(X,Y) is a row vector containing the
%   application of QUAD_OVER_LIN to each column. For N-D arrays, the
%   operation is applied to the first non-singleton dimension of X.
%
%   QUAD_OVER_LIN(X,Y,DIM) takes the sum along the dimension DIM of X.
%   A special value of DIM == 0, is accepted here, which is automatically
%   replaced with DIM == NDIMS(X) + 1. This has the effect of eliminating
%   the sum; thus QUAD_OVER_LIN( X, Y, NDIMS(X) + 1 ) = ABS( X ).^2 ./ Y.
%
%   In all cases, Y must either be a scalar or a matrix of the same size
%   as SUM(X,DIM).
%
%   Disciplined quadratic programming information:
%       QUAD_OVER_LIN is convex, nonmontonic in X, and nonincreasing in Y.
%       Thus when used with CVX expressions, X must be convex (or affine)
%       and Y must be concave (or affine).

%
% Check arguments
%

if ~isreal( y ),
    error( 'Second argument must be real.' );
elseif nargin < 3 | isempty( dim ),
    dim = cvx_default_dimension( size( x ) );
elseif ~cvx_check_dimension( dim ),
    error( 'Third argument, if supplied, must be a positive integer.' );
end

%
% Perform calculation
%

z = sum_square( x, dim );
if length( y ) ~= 1 & ~isequal( size( z ), size( y ) ),
    error( 'Input size mismatch.' );
end
temp = y <= 0;
inf_fix = any( temp );
if inf_fix,
    y( temp ) = 1;
end
z = z ./ y;
if inf_fix,
    z( temp ) = +Inf;
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
