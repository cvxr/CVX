function y = cvxobj( p )
global cvx___
if isempty( cvx___ ), 
    error( 'Internal cvx data corruption' ); 
elseif isa( p, 'cvxprob' ),
    p = index( p );
elseif ~isnumeric( p ) | length( p ) ~= 1 | p <= 0 | p > length( cvx___.problems ) + 1 | p ~= floor( p ),
    error( 'Argument must be a positive integer.' );
end
y = class( struct( 'index_', p, 'id_', cvx___.id + 1 ), 'cvxobj' );
cvx___.id = y.id_ + 1;

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
