function y = horzcat( varargin )
error( cvx_verify( varargin{:} ) );
y = cat( 2, varargin{:} );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
