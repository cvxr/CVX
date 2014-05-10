function y = cvxdual( prob, name )
if ischar( name ), name = struct( 'type', '.', 'subs', name ); end
y = class( struct( 'problem_', prob, 'name_', name, 'attached_', false, 'id_', cvx_id( cvx ) ), 'cvxdual' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
