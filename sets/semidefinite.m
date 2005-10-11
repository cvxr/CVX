function cvx_optpnt = semidefinite( n, iscplx )
error( nargchk( 1, 2, nargin ) );

%
% Check size vector
%

if ~isnumeric( n ) | length( n ) ~= 1 | n <= 0 | n ~= floor( n ),
    error( 'First argument must be a positive integer.' );
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

cvx_begin_set
   if iscplx,
       variable x( n, n ) hermitian
   else,
       variable x( n, n ) symmetric
   end
   global cvx___
   p = index( cvx_problem );
   nn = 2 : length( cvx___.problems( p ).reserved );
   cvx___.problems( p ).reserved( nn ) = 1;
   if iscplx,
       s = 'hermitian-semidefinite';
   else,
       s = 'semidefinite';
   end
   cvx___.problems( p ).cones = struct( 'type', s, 'indices', nn( : ) );
cvx_end_set

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
