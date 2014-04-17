function y = cvxobj()
persistent id
if isempty( id ), id = 0; end
id = id + 1;
y = class( struct( 'id_', id ), 'cvxobj' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
