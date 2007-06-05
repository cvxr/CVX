function cvx_begin( varargin )

% CVX_BEGIN    Starts a new cvx specification.
%
% The CVX_BEGIN command is used to mark the beginning of a cvx model.
% Following this command may be variable declarations, objective functions,
% and constraints, and a CVX_END to mark the completion of the model.
%
% If another model has already been created and is still in progress, then
% CVX_BEGIN will abort with an error.
%
% The command
%    CVX_BEGIN SDP
% can be used to mark the beginning of a semidefinite programming (SDP)
% model in CVX. This command alters the interpretation of inequality
% constraints when used with matrices, so that SDPs are easier to construct.
% Specifically, for constraints of the form
%     X >= Y    X > Y    Y < X    Y <= X
% where X and Y are matrices (i.e., not vectors or scalars), CVX will
% interpret them all as LMIs, and convert them to
%     X - Y == semidefinite(size(X,1));
% X and Y _must_ be square and identically sized---with one exception: 
% X or Y may also be the scalar number zero, so that expressions such as
% X >= 0 have the expected meaning.
%
% The command
%   CVX_BEGIN GP
% can be used to mark the beginning of a geometric programming (GP) model
% in CVX. This command alters the definition of the VARIABLE keyword so
% that it creates geometric variables by default. GP and SDP cannot be
% supplied for the same problem.
%
% The command
%    CVX_BEGIN SET
% can be used to mark the beginning of a set definition---a cvx feasibility
% problem intended to describe a set for use in other models. See the files
% in the cvx subdirectory sets/ for examples. The SET keyword can be
% combined with SDP or GP to specify sets which use SDP or GP conventions;
% for example,
%   CVX_BEGIN SET SDP
%
% The command
%     CVX_BEGIN SEPARABLE
% gives permission for CVX to solve a multiobjective problem simply by
% taking the sum of the objectives and solving the resulting single-
% objective problem. As the name implies, this produces equivalent results
% when the subproblems are separable---and when there is not a mixture of
% infeasible, unbounded, and feasible subproblems. (A later version of CVX
% will remove this latter limitation.) This keyword is useful, for example,
% when the CVX model is being used to compute a scalar function applied
% elementwise to an array. It is ignored for single-objective problems,
% feasibility problems, sets, and incomplete specifications.
%
% The command
%     CVX_BEGIN SEPARABLE
% gives permission for CVX to solve a multiobjective problem simply by
% taking the sum of the objectives and solving the resulting single-
% objective problem. As the name implies, this produces equivalent results
% when the subproblems are separable---and when there is not a mixture of
% infeasible, unbounded, and feasible subproblems. (A later version of CVX
% will remove this latter limitation.) This keyword is useful, for example,
% when the CVX model is being used to compute a scalar function applied
% elementwise to an array. It is ignored for single-objective problems,
% feasibility problems, sets, and incomplete specifications.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

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
