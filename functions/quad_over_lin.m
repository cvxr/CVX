function cvx_optval = quad_over_lin( x, y, dim )
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

sx = size( x );
if nargin < 3 | isempty( dim ),
    dim = [ find( sx ~= 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | dim < 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension.' );
elseif dim == 0,
    dim = find( sx == 1 );
    if isempty( dim ),
        dim = length( sx ) + 1;
    else,
        dim = dim( 1 );
    end
end

%
% Check sizes
%

sy = size( y );
nd = max( [ length( sx ), length( sy ), dim ] );
sx = [ sx, ones( 1, nd - length( sx ) ) ];
sy = [ sy, ones( 1, nd - length( sy ) ) ];
sz = sx;
sz( dim ) = 1;
if any( sy ~= 1 ) & any( sy ~= sz ),
    error( 'Dimensions are not compatible.' );
end

%
% Check curvature
%

if ~cvx_isaffine( x ),
    error( 'The first argument must be affine.' );
elseif ~cvx_isconcave( y ),
    error( 'The second argument must be concave.' );
end

%
% Construct the problem. Note that in order to support the case
% where y is convex, we actually have to introduce a new inequality
% y2 <= y, where y2 is a new variable. That is because set membership
% constraints require affine arguments. So this is an obscure case, then,
% where automatic curvature checking in epigraph functions breaks down.
%

cvx_begin
    variable z( sz )
    minimize z
    if cvx_isconcave( y ),
        variable y2( sy )
	    y2 <= y;
	else
	    y2 = y;
	end
    { x, y2, z } == rotated_lorentz( sx, dim, ~isreal( x ) );
cvx_end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
