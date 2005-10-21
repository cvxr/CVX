function z = newslack( prob, sz, str )
error( nargchk( 1, 3, nargin ) );
if nargin < 3, str = []; end

global cvx___
p = index( prob );
base( 1 ).type = '.';
base( 1 ).subs = 'slack_';
base( 2 ).type = '{}';
base( 2 ).subs = { eval( 'length(cvx___.problems( p ).variables.slack_)+1','1' ) };
z = newvar( prob, base, sz, str );
temp = find( any( cvx_basis( z ), 1 ) );
temp = temp( : )';
cones = cvx___.problems( p ).cones;
ncone = struct( 'type', 'nonnegative', 'indices', temp );
if isempty( cones ),
    cones = ncone;
elseif isequal( cones( 1 ).type, 'nonnegative' ),
    cones( 1 ).indices = [ cones( 1 ).indices, temp ];
else,
    cones = [ ncone, cones ];
end
cvx___.problems( p ).cones = cones;
cvx___.problems( p ).reserved( temp ) = 1;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
