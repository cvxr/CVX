function cvx_begin( varargin )

%CVX_BEGIN    Starts a new cvx specification.
%   CVX_BEGIN marks the beginning of a new cvx model. Following this command
%   may be variable declarations, objective functions, and constraints, and
%   a CVX_END to mark the completion of the model.
%
%   If another model has already been created and is still in progress, then
%   CVX_BEGIN will issue a warning and clear the previous model.
%
%   CVX_BEGIN SDP marks the beginning of a semidefinite programming (SDP) model.
%   This command alters the interpretation of inequality constraints when used
%   with matrices, so that SDPs are easier to construct. Specifically, 
%   constraints of the form
%       X >= Y    X > Y    Y < X    Y <= X
%   where X and Y are matrices (i.e., not vectors or scalars), CVX will
%   interpret them all as LMIs, and convert them to
%        X - Y == semidefinite(size(X,1));
%   X and Y _must_ be square and identically sized---with one exception: 
%   X or Y may also be the scalar number zero, so that expressions such as
%   X >= 0 have the expected meaning.
%
%   CVX_BEGIN GP marks the beginning of a geometric programming (GP) model. This
%   command alters the definition of the VARIABLE keyword to create geometric
%   variables by default GP and SDP cannot be supplied simultaneously.
%
%   CVX_BEGIN SET can be used to mark the beginning of a set definition---a cvx
%   feasibility problem intended to describe a set for use in other models. See
%   the files in the cvx subdirectory sets/ for examples. The SET keyword can be
%   combined with SDP or GP to specify sets which use SDP or GP conventions;
%   for example, CVX_BEGIN SET SDP
%
%   CVX_BEGIN SEPARABLE gives permission for CVX to solve a multiobjective
%   problem simply by taking the sum of the objectives and solving the resulting
%   single-objective problem. As the name implies, this produces equivalent 
%   results only when the subproblems are separable. Behavior is undefined when
%   one or more of the subproblems is infeasible or unbounded. The keyword is
%   ignored for sets and incomplete specifications.

if ~iscellstr( varargin ),
    error( 'Arguments must be strings.' );
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    if ~isempty( cvx_problem.objective ) | ~isempty( cvx_problem.variables ) | ~isempty( cvx_problem.duals ) | nnz( cvx_problem.t_variable ) > 1,
        warning( sprintf( 'A cvx problem already existed in this scope.\n   It is being overwritten.' ) );
    end
    evalin( 'caller', 'pop( cvx_problem, ''clear'' );' );
    cvx_problem = [];
%   error( sprintf( 'A cvx problem already exists in this scope.\n(To clear it and start a new one, use the command ''cvx_clear''.)' ) );
end
cvx_setpath( 1 );
global cvx___
if isempty( cvx___.problems ) & cvx___.profile,
    profile resume
end
cvx_create_problem( varargin{:} );
assignin( 'caller', 'cvx_problem', cvx_problem );

% Copyright 2007 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.