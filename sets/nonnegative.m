function cvx_optpnt = nonnegative( sx )
error( nargchk( 1, 1, nargin ) );

[ temp, sx ] = cvx_check_dimlist( sx, false );
if ~temp,
    error( 'Argument must be a non-empty dimension vector.' );
end
    
cvx_begin_set
   variables x( sx )
   global cvx___
   p = index( cvx_problem );
   temp = [ 2 : length( cvx___.problems( p ).reserved ) ]';
   cvx___.problems( p ).reserved( temp ) = 1;
   cvx___.problems( p ).cones = struct( 'type', 'nonnnegative', 'indices', temp );
cvx_end_set

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
