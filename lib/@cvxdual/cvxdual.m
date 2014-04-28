function y = cvxdual( prob, name )

global cvx___
id = cvx___.id + 1;
cvx___.id = id;
y = class( struct( 'problem_', prob, 'name_', name, 'attached_', false, 'id_', id ), 'cvxdual' );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
