function minimize( varargin )

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isempty( cvx_problem ) | ~isa( cvx_problem, 'cvxprob' ),
    error( 'A cvx problem does not exist in this scope.' );
end
if ~isempty( cvx_problem.direction ),
    if isequal( cvx_problem.direction, 'find' ),
        error( 'Objective functions cannot be added to sets.' );
    else,
        error( 'An objective function has already been supplied.' );
    end
end

if nargin < 1,
    error( 'Objective expression missing.' );
elseif iscellstr( varargin ),
    varargin = { evalin( 'caller', sprintf( '%s ', varargin{:} ) ) };
end

if ~cvx_isvalid( varargin ),
    error( 'Input argument(s) must be valid cvx objective expressions.' );
elseif length( varargin ) > 1 | isa( varargin{1}, 'cell' ) | isa( varargin{1}, 'struct' ),
    error( 'Objective functions cannot be composite.' );
end

if ~cvx_isconvex( varargin{1} ),
    error( sprintf( 'Dispciplined convex programming error:\n   Objective function in a minimization must be convex.' ) );
end

newobj( cvx_problem, 'minimize', varargin{1} );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
