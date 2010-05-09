function sout = cvx_solver( flag )

%CVX_SOLVER    CVX solver selection.
%   CVX_SOLVER SDPT3  selects SDPT3  as the solver CVX uses to solve models.
%   CVX_SOLVER SeDuMi selects SeDuMi as the solver CVX uses to solve models.
%   The solver string is case-insensitive; e.g., 'sdpt3' is fine.
%
%   The default solver is SeDuMi. For most models, this will be the best
%   choice. Nevertheless, for your problems you may find that SDPT3 is
%   more accurate or reliable, hence the need for this command.
%
%   If CVX_SOLVER is called within a model---that is, between the statements
%   CVX_BEGIN and CVX_END---then the new solver selection applies only to that
%   particular model. In this mode, the change takes place for *both*
%   successive approximation models and normal models.
%
%   If called outside of a model, then the change applies to
%   all subsequent models; that is, it modifies the default solver.
%
%   On exit, CVX_SOLVER returns the *previous* solver. This can be used to 
%   restore that solver selection at a later time, in a manner similar to 
%   CVX_PRECISION.
%
%   CVX_SOLVER, with no arguments, returns the current solver.

global cvx___
cvx_global
if nargin > 0,
    if isempty( flag ),
        flag = 'sedumi';
    elseif ~ischar( flag ) || size( flag, 1 ) ~= 1,
        error( 'Argument must be a string.' );
    elseif ~isfield( cvx___.path.solvers, lower( flag )  ),
        error( 'Unknown, unusable, or missing solver: %s', flag );
    end
end
if isempty( cvx___.problems ),
    s = cvx___.solver;
    if nargin > 0,
        cvx___.solver = flag;
    end
else
    s = cvx___.problems(end).solver;
    if nargin > 0,
        if ~strcmpi( s, flag ) && ~isa( evalin( 'caller', 'cvx_problem', '[]' ), 'cvxprob' ),
            warning( 'CVX:Solver', 'The global CVX solver selection cannot be changed while a model is being constructed.' );
        else
            cvx___.problems(end).solver = flag;
        end
    end
end
if nargin > 0 && cvx___.path.hold,
    cvx_setspath( flag );
end
if nargin == 0 || nargout > 0,
    sout = s;
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
