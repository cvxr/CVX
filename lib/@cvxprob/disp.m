function disp( prob, varargin )
cvx_display( cvx_validate( prob.index_, prob.id_ ), varargin{:} );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
