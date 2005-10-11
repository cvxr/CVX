function cvx_optval = sum_largest( x, k, dim )
error( nargchk( 2, 3, nargin ) );

%SUM_LARGEST   Sum of the largest k elements of an array.
%
%   For a real vector X and an integer k between 1 and length(X) inclusive,
%   y = SUM_LARGEST(X,k) is the sum of the k largest elements of X; e.g.,
%       temp = sort( x )
%       y = sum( temp( 1 : k ) )
%   If k=1, then SUM_LARGEST(X,k) is equivalent to MAX(X); if k=length(X),
%   then SUM_LARGEST(X,k) is equivalent to SUM(X).
%
%   Both X and k must be real, and k must be a scalar. But k is not, in
%   fact, constrained to be an integer between 1 and length(X); the
%   function is extended continuously and logically to all real k. For
%   example, if k <= 0, then SUM_LARGEST(X,k)=0. If k > length(X), then
%   SUM_LARGEST(X,k)=SUM(X). Non-integer values of k interpolate linearly
%   between their integral neighbors.
%
%   For matrices, SUM_LARGEST(X,k) is a row vector containing the
%   application of SUM_LARGEST to each column. For N-D arrays, the 
%   SUM_LARGEST operation is applied to the first non-singleton dimension
%   of X.
%
%   SUM_LARGEST(X,k,DIM) performs the operation along dimension DIM of X.
%
%   Disciplined convex programming information:
%       SUM_LARGEST(X,...) is convex and nondecreasing in X. Thus, when
%       used in CVX expressions, X must be convex (or affine). k and DIM
%       must both be constant.

sx = size( x );
if nargin < 3 | isempty( dim ),
    dim = [ find( sx > 1 ), 1 ];
    dim = dim( 1 );
elseif ~isnumeric( dim ) | ~isreal( dim ) | dim <= 0 | dim ~= floor( dim ),
    error( 'Second argument must be a dimension.' );
elseif ~isnunmeric( k ) | ~isreal( k ),
    error( 'Third argument must be a real scalar.' );
elseif ~isreal( x ),
    error( 'First argument must be real.' );
end

%
% Quick exit cases
%

nd = max( dim, length( sx ) );
sx = [ sx, ones( 1, dim - nd ) ];
sy = sx;
sy( dim ) = 1;
if k <= 0,
    
    cvx_optval = zeros( sy );
    
elseif k <= 1,
    
    cvx_optval = k * max( x, [], dim );
    
elseif k >= sx( dim ),
    
    cvx_optval = sum( x, dim );
    
else,
    
    nd = length( sx );
    ndxs = cell( 1, nd );
    [ ndxs{:} ] = deal( ':' );
    ndxs{ dim } = ones( 1, sx( dim ) );
    
    cvx_begin
        variables xp( sx ) yp( sy ) z( sy )
        minimize z
        subject to
            xp >= 0;
            k * yp == sum( xp, dim ) - z;
            xp >= cvx_subref( yp, ndxs{:} ) + x;
    cvx_end
    
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
