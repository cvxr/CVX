function cvx_optval = geomean( x, dim )
error( nargchk( 1, 2, nargin ) );

% GEOMEAN   Geometric mean.
%     For vectors, GEOMEAN(X) is the geometric mean of the elements of X.
%     If any of the elements of X are negative, then Y=-Inf. Otherwise,
%     Y = PROD(X)^(1/LENGTH(X)). All elements must be real.
%
%     For matrices, GEOMEAN(X) is a row vector containing the geometric
%     means of the columns. For N-D arrays, GEOMEAN(X) is an array of the
%     geometric means taken along the first non-singleton dimension of X.
%
%     GEOMEAN(X,DIM) takes the geometric mean along the dimension DIM of X.
%
%     For matrices, Y=GEOMEAN(X) is a vector of the geometric means of the
%     columns of X. For N-D arrays, the geometric means are taken along the
%     first non-singleton dimension of X.
%
%     Disciplined convex programming information:
%         GEOMEAN is convex and nondecreasing; therefore, when used in CVX
%         specifications, its argument must be concave.

%
% Check types
%

if ~isreal( x ),
    error( 'First argument must be real.' );
elseif ~cvx_isconcave( x ),
    error( sprintf( 'Disciplined convex programming error: \n    GEOMEAN is concave and nondecreasing; its argument must therefore be concave.' ) );
elseif nargin < 2,
    dim = cvx_default_dimension( size( x ) );
elseif ~cvx_check_dimension( dim, false ),
    error( 'Second argument must be a positive integer.' );
end

%
% Determine sizes
%

sx = size( x );
sx = [ sx, ones( 1, dim - length( sx ) ) ];
nx = sx( dim );

%
% Empty arrays
%

if any( sx == 0 ),
    sx( dim ) = 1;
    cvx_optval = ones( sx );
    return
end

%
% A single element: -Inf if x < 0, x if x >= 0
%

if nx == 1,
    cvx_begin
        variable z( sz )
        maximize z;
        z <= x;
        0 <= x;
    cvx_end
    return
end

%
% Permute the matrix, if needed, so the geometric mean can be taken
% along the first dimension.
%

if dim > 1 & any( sx( 1 : dim - 1 ) > 1 ),
    perm = [ dim, 1 : dim - 1, dim + 1 : length( sx ) ];
    x = permute( x, perm );
    sx = sx( perm );
    dim = 1;
else,
    perm = [];
end

%
% Construct the problem.
% --- For n == 2, we use the following equivalency
%     sqrt(x*y)>=z, x>=0, y>=0 <--> z^2/x <= y
% --- For n == 2^k, we can recursively apply this log(k,2) times.
% --- For other n, note that
%     x(1)*x(2)*...*x(n) <= y^n
%         <----> x(1)*x(2)*...*x(n)*y^m <= y^(n+m)
%     so by adding extra y's to the left-hand side we can use the same
%     recursion for lengths that are not powers of two.
%

nv = prod( sx ) / nx;
x = reshape( x, nx, nv );
cvx_begin
    tt = ~cvx_isaffine( x, true );
    if any( tt ),
        temp = x( tt );
        variable x2( size( temp ) );
        x( tt ) = x2;
        temp >= x2;
    end
    variable y( 1, nv )
    maximize( y )
    while nx >= 2,
        n2 = ceil( nx * 0.5 );
        if n2 * 2 ~= nx, x = [ x ; y ]; end
        cone = rotated_lorentz( [ n2, nv ], 0 );
        cone.y == x( 1 : 2 : end, : );
        cone.z == x( 2 : 2 : end, : );
        x = cone.x;
        nx = n2;
    end
    y == x;
cvx_end

%
% Reverse the reshaping and permutation steps
%

sx( dim ) = 1;
cvx_optval = reshape( cvx_optval, sx );
if ~isempty( perm ),
    cvx_optval = ipermute( cvx_optval, perm );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
