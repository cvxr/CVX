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
%    CVX_BEGIN SET
% can be used to mark the beginning of a set definition---a cvx feasibility
% problem intended to describe a set for use in other models. See the files
% in the cvx subdirectory sets/ for examples.
%
% The commands
%    CVX_BEGIN SET SDP
%    CVX_BEGIN SDP SET
% are also possible, and combine the two functions. That is, they mark the
% beginning of a set definition which uses the SDP constraint conventions
% for its description.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

if ~iscellstr( varargin ),
    error( 'Arguments, if supplied, must be strings.' );
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
cvx_setpath( 1 );
cvx_create_problem( varargin{:} );
assignin( 'caller', 'cvx_problem', cvx_problem );
