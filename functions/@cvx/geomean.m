function y = geomean( x, dim, w, ismap )
error( nargchk( 1, 4, nargin ) );

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
%     GEOMEAN(X,DIM,W), where W is a vector of nonnegative integers,
%     computes a weighted geometric mean Y = PROD(X.^W)^(1/SUM(W)). In 
%     effect, each element X(i) is replicated W(i) times in the geometric
%     mean---but this method is far more efficient.
%
%     Disciplined convex programming information:
%         GEOMEAN is convex and nondecreasing; therefore, when used in CVX
%         specifications, its argument must be concave.

%
% Basic argument check
%

sx = size( x );
if nargin < 2 | isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim ),
    error( 'Second argument must be a positive integer.' );
end

%
% Determine sizes, quick exit for empty arrays
%

sx = [ sx, ones( 1, dim - length( sx ) ) ];
nx = sx( dim );
sy = sx;
sy( dim ) = 1;
if any( sx == 0 ),
    y = ones( sy );
    return
end

%
% Third and fourth argument check
%

map = [];
if nargin == 4,
    map = w;
elseif nargin < 3 | isempty( w ),
    w = [];
elseif numel( w ) ~= length( w ) | ~isnumeric( w ) | ~isreal( w ) | any( w < 0 ) | any( w ~= floor( w ) ),
    error( 'Third argument must be a vector of nonnegative integers.' );
elseif length( w ) ~= nx,
    error( sprintf( 'Third argument must be a vector of length %d', nx ) );
elseif ~any( w ),
    y = ones( sy );
    return
end

%
% Type check
%

persistent remap_1 remap_2 remap_3 remap_4
if isempty( remap_4 ),
    % Constant (postive or negative)
    remap_1 = cvx_remap( 'real' );
    remap_2 = cvx_remap( 'concave' );
    remap_3 = cvx_remap( 'log-convex' );
    remap_4 = cvx_remap( 'log-concave' );
end
vx = cvx_reshape( cvx_classify( x ), sx );
t1 = all( reshape( remap_1( vx ), sx ), dim );
t2 = all( reshape( remap_2( vx ), sx ), dim );
t3 = all( reshape( remap_3( vx ), sx ), dim ) | ...
     all( reshape( remap_4( vx ), sx ), dim );
% Valid combinations with zero or negative entries can be treated as constants
t1 = t1 | ( ( t2 | t3 ) & any( vx == 1 | vx == 9, dim ) );
ta = t1 + ( 2 * t2 + 3 * t3 ) .* ~t1;
nu = unique( ta( : ) );
nk = length( nu );

%
% Permute and reshape, if needed
%

perm = [];
if nk > 1 | ( any( nu > 1 ) & nx > 1 ),
    if dim > 1 & any( sx( 1 : dim - 1 ) > 1 ),
        perm = [ dim, 1 : dim - 1, dim + 1 : length( sx ) ];
        x   = permute( x,  perm );
        sx  = sx( perm );
        sy  = sy( perm );
        ta  = permute( ta, perm );
        dim = 1;
    else
        perm = [];
    end
    nv = prod( sy );
    x  = reshape( x, nx, nv );
    ta = reshape( ta, 1, nv );
end

%
% Perform the computations
%

if nk > 1,
    y = cvx( [ 1, nv ], [] );
end
for k = 1 : nk,

    if nk == 1,
        xt = x;
        sz = sy;
    else
        tt = ta == nu( k );
        xt = cvx_subsref( x, ':', tt );
        nv = nnz( tt );
        sz = [ 1, nv ];
    end

    switch nu( k ),
        case 0,
            error( sprintf( 'Disciplined convex programming error:\n   Invalid computation: geomean( {%s} )', cvx_class( xt, true, true ) ) );
        case 1,
            yt = geomean( cvx_constant( xt ), dim, w );
        case 2,
            if nx == 1,
                cvx_begin
                    variable z( sz )
                    z <= xt;
                    z <= 0;
                cvx_end
                yt = cvx_optval;
            else
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
                % --- For integer-weighted geometric means, we effectively
                %     replicate each x(i) w(i) times. However, because
                %     sqrt(x(i)*x(i)) = x(i), we can prune away most of
                %     the duplicated values very cheaply. The number of
                %     non-trivial appearances of x(i) will be reduced from
                %     w(i) to at most ceil(log2(w(i))), or more precisely
                %     the number of 1's in a binary expansion of w(i).
                %
                if isempty( map ),
                    if isempty( w ),
                        w = ones( 1, nx );
                    end
                    map  = cvx_geomean_map( w );
                end
                nm   = size(map,2);
                msiz = [ 1, nm, nv ];
                cvx_begin
                    variable xw( nm, nv );
                    xt = [ cvx_accept_concave( xt ) ; xw ];
                    cone = rotated_lorentz( msiz, 0 );
                    cone.y == reshape( xt(map(1,:),:), msiz ); 
                    cone.z == reshape( xt(map(2,:),:), msiz );
                    cone.x == reshape( xt(map(3,:),:), msiz );
                    maximize( xw( end, : ) );
                cvx_end
                yt = cvx_optval;
            end
        case 3,
            if nx == 1,
                yt = xt;
            elseif isempty( w ),
                yt = exp( sum( log( xt ), 1 ) * ( 1 / nx ) );
            else
                yt = exp( ( w / sum( w ) ) * log( xt ) );
            end
        otherwise,
            error( 'Shouldn''t be here.' );
    end

    if nk == 1,
        y = yt;
    else
        y = cvx_subsasgn( y, tt, yt );
    end

end

%
% Reverse the reshaping and permutation steps
%

y = reshape( y, sy );
if ~isempty( perm ),
    y = ipermute( y, perm );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
