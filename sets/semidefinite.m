function cvx_optpnt = semidefinite( n, iscplx )

%SEMIDEFINITE   Real symmetric positive semidefinite matrices.
%    SEMIDEFINITE(N), where N is an integer, creates a symmetric matrix
%    variable of size [N,N] and constrains it to be positive semidefinite.
%    Therefore, given the declaration
%       variable x(n,n) symmetric
%    the constraint
%       x == semidefinite(n)
%    is equivalent to
%       lambda_min(x) >= 0;
%    In fact, lambda_min is implemented in CVX using SEMIDEFINITE for
%    real matrices.
%
%    SEMIDEFINITE(SX), where SX is a valid size vector, creates an array
%    variable of size SX and constrains each subarray along the leading two
%    dimensions to be positive semidefinite. SX(1) and SX(2) must be equal.
%    Therefore, given the declaration
%       variable x(sx) symmetric
%    the constraint
%       x == semidefinite(sx)
%    is equivalent to
%       for k = 1:prod(sx(3:end)),
%          lambda_min(x(:,:,k)) >= 0;
%       end
%
%    SEMIDEFINITE(N,CPLX) and SEMIDEFINITE(SX,CPLX) create real semidefinite
%    sets if CPLX is FALSE, and complex Hermitian semidefinite sets if CPLX
%    is TRUE. The latter case is equivalent to calling the function
%    HERMITIAN_SEMIDEFINITE.
%
%   Disciplined convex programming information:
%       SEMIDEFINITE is a cvx set specification. See the user guide for
%       details on how to use sets.

%
% Check size vector
%

error( nargchk( 1, 2, nargin ) );
if ~isnumeric( n ) | isempty( n ) | any( n < 0 ) | any( n ~= floor( n ) ),
    error( 'First argument must be a positive integer or a valid size vector.' );
elseif length( n ) > 1 & n( 1 ) ~= n( 2 ),
    error( 'If a size vector is supplied, the first two dimensions must be identical.' );
end

%
% Check iscplx flag
%

if nargin < 2,
    iscplx = false;
elseif ( ~isnumeric( iscplx ) & ~islogical( iscplx ) ) | length( iscplx ) ~= 1,
    error( 'Second argument must be a numeric or logical scalar.' );
end

%
% Construct the cone
%

if length( n ) == 1,
    sz = [ n, n ];
    nv = 1;
else
    sz = n;
    n = sz( 1 );
    nv = prod( sz( 3 : end ) );
end
if iscplx,
    ntri = n * n;
else
    ntri = n * ( n + 1 ) / 2;
end
if ntri == 0 | nv == 0,
    cvx_optpnt = zeros( sz );
    return
end
cvx_begin_set
   if iscplx,
       variable x( sz ) hermitian
   else
       variable x( sz ) symmetric
   end
   global cvx___
   p = index( cvx_problem );
   if iscplx,
       s = 'hermitian-semidefinite';
   else
       s = 'semidefinite';
   end
   tx = reshape( find( any( cvx_basis( x ), 2 ) ), ntri, nv );
   newnonl( cvx_problem, s, tx );
cvx_end_set

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
