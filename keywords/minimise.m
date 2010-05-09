function minimise( varargin )

%MINIMISE Specifiies a convex (or affine) objective to be maximized.

cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isempty( cvx_problem ) || ~isa( cvx_problem, 'cvxprob' ),
    error( 'A cvx problem does not exist in this scope.' );
end
if ~isempty( cvx_problem.direction ),
    if isequal( cvx_problem.direction, 'find' ),
        error( 'Objective functions cannot be added to sets.' );
    else
        error( 'An objective function has already been supplied.' );
    end
end

if nargin < 1,
    error( 'Objective expression missing.' );
elseif iscellstr( varargin ),
    arg = evalin( 'caller', sprintf( '%s ', varargin{:} ) );
elseif nargin > 1,
    error( 'Too many input arguments.' );
else
    arg = varargin{1};
end

if ~isa( arg, 'cvx' ) && ~isa( arg, 'double' ) && ~isa( arg, 'sparse' ),
    error( 'Cannot accept an objective of type ''%s''.', class( arg ) );
end
persistent remap
if isempty( remap ),
    remap = cvx_remap( 'convex', 'log-convex' );
end
vx = remap( cvx_classify( arg ) );
if ~all( vx ),
    error( 'Disciplined convex programming error:\n   Cannot minimise a(n) %s expression.', cvx_class(arg(vx==0),false,true) );
end

newobj( cvx_problem, 'minimize', arg );

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
