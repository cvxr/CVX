function cvx_optval = max( x, y, dim )
error( nargchk( 1, 3, nargin ) );

% MAX    Largest component (augmented with CVX support).
%    For vectors, MAX(X) is the largest element in X. For matrices,
%    MAX(X) is a row vector containing the maximum element from each
%    column. For N-D arrays, MAX(X) operates along the first
%    non-singleton dimension.
% 
%    MAX(X,[],DIM) operates along the dimension DIM.
%
%    MAX(X,Y) returns an array the same size as X and Y with the
%    largest elements taken from X or Y.  Either one can be a scalar.
%
%    For CVX expresisons, the inputs X and/or Y must be real. For numeric
%    expressions, they can be complex; the magnitude is used, and the angle
%    is ignored. NaN's are ignored when computing the maximum.
%
%    Disciplined convex programming information:
%        MAX is convex and nondecreasing in its first two arguments. Thus
%        when used in CVX expressions, those inputs must be convex (or
%        affine). The third expression must be constant.

if nargin == 2,

    %
    % max( X, Y )
    %

    if ~cvx_isconvex( x ) | ~cvx_isconvex( y ),
        error( sprintf( 'Disciplined convex programming error:\n   max() is convex and nondecreasing, so its arguments must be convex.' ) );
    end

    sx = size( x );
    sy = size( y );
    if all( sx == 1 ),
        sz = sy;
    elseif all( sy == 1 ),
        sz = sx;
    elseif ~isequal( sx, sy ),
        error( 'Array dimensions must match.' );
    end

    %
    % The cvx multi-objective problem
    %

    cvx_begin
        variable z( sz );
        minimize z
        subject to
            x <= z;
            y <= z;
    cvx_end
    
else,
    
    %
    % max( X, [], dim )
    %

    if ~cvx_isconvex( x ),
        error( sprintf( 'Disciplined convex programming error:\n   max() is convex and nondecreasing, so its argument must be convex.' ) );
    elseif nargin > 1 & ~isempty( y ),
        error( 'max with two matrices to compare and a working dimension is not supported.' );
    end
        
    sx = size( x );
    if nargin == 1,
        dim = [ find( sx > 1 ), 1 ];
        dim = dim( 1 );
    elseif ~isnumeric( dim ) | dim < 0 | dim ~= floor( dim ),
        error( 'Dimension argument must be a positive integer.' );
    end

    %
    % Trivial exit if the dimension has no width
    %

    if dim > length( sx ) | sx( dim ) == 1,
        z = x;
        return
    end

    %
    % Construct the (possibly multiobjective) problem
    %

    cvx_begin
        variable z( sx( 1 : dim - 1 ), 1, sx( dim + 1 : end ) );
            z2 = reshape( z, [ prod( sx( 1 : dim - 1 ) ), 1, prod( sx( dim + 1 : end ) ) ] );
            z2 = cvx_subref( z2, ':', ones( 1, sx( dim ) ), ':' );
            z2 = reshape( z2, sx );
        minimize z
        subject to
            x <= z2;
    cvx_end
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
