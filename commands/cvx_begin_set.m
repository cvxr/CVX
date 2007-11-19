function cvx_begin_set( varargin )

% CVX_BEGIN_SET    Starts a new cvx CVX specification.
%
% The command
%    CVX_BEGIN_SET
% can be used to mark the beginning of a set definition---a cvx feasibility
% problem intended to describe a set for use in other models. See the files
% in the cvx subdirectory sets/ for examples. It has actually been
% deprecated in favor of
%    CVX_BEGIN SET
% (two separate words). However, it will continue to be available for back-
% compatability reasons, so you are free to use either.

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
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
    cvx_pop( cvx_problem );
end
cvx_global
cvx_setpath( 1 );
if cvx___.profile & isempty( cvx___.problems ),
    profile resume
end
cvx_create_problem( 'set', varargin{:} );
assignin( 'caller', 'cvx_problem', cvx_problem );
