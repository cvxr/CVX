function varargout = in( varargin )
error( cvx_verify( varargin{:} ) );

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if ~isa( cvx_problem, 'cvxprob' ),
    error( 'A valid cvx problem must be created first.' );
else,
    eq( varargin{:} );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
