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
%   A special value of DIM == 0, is accepted here, has the effect of
%   eliminating the sum; thus QUAD_OVER_LIN(X,Y,0) = ABS(X).^2 ./ Y.
%
%   In all cases, Y must be real and have a size compatible with elementwise
%   division with SUM(X,DIM); and in particular, if DIM ~= 0, then
%   size(Y,DIM) == 1.
%
%   Disciplined quadratic programming information:
%       QUAD_OVER_LIN is convex, nonmontonic in X, and nonincreasing in Y.
%       Thus when used with CVX expressions, X must be convex (or affine)
%       and Y must be concave (or affine).

sx = size( x );
sy = size( y );
nx = length( sx );
ny = length( sy );
nd = max( nx, ny );
if nargin < 3 | isempty( dim ),
    dim = [ find( sx ~= 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | dim < 0 | dim ~= floor( dim ),
    error( 'Third argument must be a dimension.' );
elseif dim == 0,
    dim = nd + 1;
    nd = dim;
end

%
% Check sizes
%

sx = [ sx, ones( 1, nd - nx ) ];
sy = [ sy, ones( 1, nd - ny ) ];
sz = sx;
sz( dim ) = 1;
if sy( dim ) ~= 1 & any( sx ~= 1 ),
    error( 'Dimensions are not compatible.' );
end
if all( sz == sy ) | all( sy == 1 ),
    need_contraction = false;
elseif all( sz == 1 ),
    sz = sy;
    need_contraction = sx( dim ) ~= 1;
    sx = sz;
    dim = 0;
else,
    error( 'Dimensions are not compatible.' );
end

%
% Check curvature
%

if ~cvx_isaffine( x ),
    error( 'The first argument must be affine.' );
elseif ~isreal( y ),
    error( 'The second argument must be real.' );
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

if any( sz == 0 ),
    cvx_optval = zeros( sz );
elseif any( sx == 0 ),
    cvx_begin
        y == nonnegative( sz );
    cvx_end
else,
    cvx_begin
        variable z( sz )
        minimize z
        y = cvx_accept_concave( y );
        if need_contraction,
            x = cvx_accept_convex( norms( x, 2, dim ) );
        end
        { x, y, z } == rotated_lorentz( sx, dim, ~isreal( x ) );
    cvx_end
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
