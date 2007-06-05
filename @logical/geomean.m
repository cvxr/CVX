function y = geomean( x, dim )
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
% Check arguments
%

if ~isreal( x ), 
    error( 'First argument must be real.' ); 
elseif nargin < 2,
    dim = cvx_default_dimension( size( x ) );
elseif ~cvx_check_dimension( dim ),
    error( 'Second argument must be a positive integer.' );
end

y = exp( sum( log( max( x, realmin ) ), dim ) * ( 1 / size( x, dim ) ) );
xmin = min( x, [], dim );
y( xmin <  0 ) = -Inf;
y( xmin == 0 ) = 0;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
