function y = size( x, varargin )
y = size( cvxaff( x ), varargin{:} );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
