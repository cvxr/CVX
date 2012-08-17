function [ sout, slist ] = cvx_solver( sname )

%CVX_SOLVER    CVX solver selection.
%   CVX_SOLVER <solver_name> or CVX_SOLVER('<solver_name>')
%   selects the named solver the CVX uses to solve models. The solver name
%   is case-insensitive; so, for example, both 'SeDuMi' and 'sedumi' will
%   select the same solver.
%
%   When CVX is first installed, the solver SDPT3 is selected as the
%   default. For most problems, this will be a good choice; nevertheless,
%   no solver is perfect, so if you encounter issues you may wish to
%   experiment with other solvers.
%
%   There are two ways to use the CVX_SOLVER command. If you use it within
%   a model---that is, between the statements CVX_BEGIN and CVX_END---then
%   the new solver selection will apply only to that particular model. For
%   instance, if the default solver is SDPT3, then the following structure
%   will solve a single model using SeDuMi instead:
%       cvx_begin
%           cvx_solver sedumi
%           variables ...
%           ...
%       cvx_end
%   On the other hand, if CVX_SOLVER is called *outside* of a model, then
%   the change will apply for all subsequent models, or until you call 
%   CVX_SOLVER once again.
%
%   [ SOLVER, SOLVER_LIST ] = CVX_SOLVER returns the name of the current
%   solver, and a cell array containing the names of all available choices.
%
%   Calling CVX_SOLVER with no input or output arguments produces a listing
%   of the solvers that CVX currently recognized, and an indication of the
%   current solver selection and/or the default.

global cvx___
cvx_global
if nargin,
    if isempty( sname ),
        sname = 'default';
    elseif ~ischar( sname ) || size( sname, 1 ) ~= 1,
        error( 'Argument must be a string.' );
    end
    try
        snumber = cvx___.solvers.map.(lower(sname));
    catch %#ok
        error( 'Unknown, unusable, or missing solver: %s', sname );
    end
    if ~isempty( cvx___.solvers.list(snumber).error ),
        error( 'Solver unusable due to prior errors: %s', sname );
    end
    if isempty( cvx___.problems ),
        cvx___.solvers.selected = snumber;
        if cvx___.solvers.active,
            cvx_setspath( snumber );
        end
    elseif ~isa( evalin( 'caller', 'cvx_problem', '[]' ), 'cvxprob' ),
        error( 'The global CVX solver selection cannot be changed while a model is being constructed.' );
    else
        cvx___.problems(end).solver.index = snumber;
    end
elseif nargout == 0,
    statvec = [ 0, cvx___.solvers.map.default, cvx___.solvers.active ];
    statstr = { 'selected', 'default', 'active' };
    if ~isempty( cvx___.problems ),
        statvec(1) = cvx___.problems(end).solver.index;
    else
        statvec(1) = cvx___.solvers.selected;
    end
    fprintf( '\n' );
    fprintf( 'Available solvers:\n' );
    fprintf( '------------------\n' );
    for k = 1 : length(cvx___.solvers.names),
        fprintf( ' %s', cvx___.solvers.names{k} );
        nstat = statstr(k==statvec);
        if ~isempty(nstat),
            nstat = sprintf( '%s,', nstat{:} );
            fprintf( ' (%s)\n', nstat(1:end-1) );
        else
            fprintf( '\n' );
        end
    end
    fprintf( '\n' );
end
if nargout > 0,
    sout = cvx___.solvers.list(cvx___.solvers.selected).name;
    slist = cvx___.solvers.names;
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
