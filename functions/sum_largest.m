function y = sum_largest( varargin )

%SUM_LARGEST Sum of the largest k values of a vector.
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

%SUM_LARGEST   Internal cvx version.

persistent P
if isempty( P ),
    P.map = cvx_remap( { 'real' ; 'convex' } );
    P.funcs = { @sum_largest_cnst, @sum_largest_cvx };
    P.constant = 1;
    P.zero = 0;
    P.reduce = true;
    P.reverse = false;
    P.name = 'sum_largest';
    P.dimarg = 3;
end
[ sx, x, k, dim ] = cvx_get_dimension( varargin, 3 );
if nargin < 2,
    cvx_throw( 'Not enough arguments.' );
elseif ~isnumeric( k ) || numel(k) ~= 1 || ~isreal( k ),
    cvx_throw( 'Second argument must be real.' );
elseif k <= 0,
    sx( dim ) = 1;
    y = zeros( sx );
    if isa( x, 'cvx' ), y = cvx( x ); end
elseif k >= sx( dim ),
    y = sum( x, dim );
elseif k <= 1,
    y = k * max( x, [], dim );
else
    y = cvx_reduce_op( P, x, k, dim );
end

function y = sum_largest_cnst( x, k )
y = sort( x, 1, 'descend' );
y = sum( y(1:floor(k),:), 1 ) + (k-floor(k)) * y(ceil(k),:);

function z = sum_largest_cvx( x, k )
persistent nneg npos
if isempty( nneg ),
    nneg = cvx_remap( 'nonnegative', 'p_nonconst' );
    npos = cvx_remap( 'nonpositive', 'n_nonconst' );
end
s = size(x);
vx = cvx_classify( x );
nn = any(reshape(nneg(vx),s),1);
np = all(reshape(npos(vx),s),1);
z = []; xp = []; yp = [];
cvx_begin
    epigraph variable z( 1, s(2) )
    variables xp( s ) yp( 1, s(2) )
    z == sum( xp, 1 ) - k * yp; %#ok
    xp >= repmat( yp, [s(1),1] ) + x; %#ok
    xp >= 0; %#ok
    cvx_setnneg( cvx_fastref( z, nn ) );
    cvx_setnpos( cvx_fastref( z, np ) );
cvx_end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

