function sout = cvx_solver( flag )

%CVX_SOLVER    CVX solver selection.
%   CVX_SOLVER SDPT3  selects SDPT3  as the solver CVX uses to solve models.
%   CVX_SOLVER SeDuMi selects SeDuMi as the solver CVX uses to solve models.
%   The solver string is case-insensitive; e.g., 'sdpt3' is fine.
%
%   The default solver is SDPT3. SeDuMi is often faster; unfortunately, it has
%   significant difficulties with problems involving second-order cones, which
%   are required for quadratics, absolute values, norms, and powers.
%
%   If CVX_SOLVER is called within a model---that is, between the statements
%   CVX_BEGIN and CVX_END---then the new solver selection applies only to that
%   particular model. If called outside of a model, then the change applies to
%   all subsequent models; that is, it modifies the default solver.
%
%   On exit, CVX_SOLVER returns the *previous* solver. This can be used to 
%   restore that solver selection at a later time, in a manner similar to 
%   CVX_PRECISION.
%
%   CVX_SOLVER, with no arguments, returns the current solver.

global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
if nargin > 0,
    if isempty( flag ),
        flag = 'SDPT3';
    elseif ~ischar( flag ) | size( flag, 1 ) ~= 1,
        error( 'Argument must be a string.' );
    elseif ~isfield( cvx___.path.solvers, lower( flag )  ),
        error( sprintf( 'Unknown, unusable, or missing solver: %s', flag ) );
    end
end
cvx_problem = evalin( 'caller', 'cvx_problem', '[]' );
if isa( cvx_problem, 'cvxprob' ),
    s = cvx_problem.solver;
    if nargin > 0,
        cvx___.problems(index(cvx_problem)).solver = flag;
    end
else
    s = cvx___.solver;
    if nargin > 0,
        cvx___.solver = flag;
    end
end
if nargin == 0 | nargout > 0,
    sout = s;
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
