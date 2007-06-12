function subject( varargin )

if nargin >= 1 & ischar( varargin{1} ) & strcmpi( varargin{1}, 'to' ),
    varargin(1) = [];
end

if length( varargin ) > 0 & iscellstr( varargin ),
    evalin( 'caller',  sprintf( '%s ', varargin{:} ) );
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
