function y = cvxdual( prob, name )

y = class( struct( 'problem_', prob, 'name_', name, 'attached_', false ), 'cvxdual', cvxobj );

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
