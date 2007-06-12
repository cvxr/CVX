function cvx_optpnt = semidefinite( n, iscplx )
error( nargchk( 1, 2, nargin ) );

%
% Check size vector
%

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
