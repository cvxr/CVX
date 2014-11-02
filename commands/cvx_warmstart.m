function res = cvx_warmstart( flag )

%CVX_WARMSTART    Controls CVX hot start capability.
%   The CVX_WARMSTART command controls CVX's warm start functionality.
%
%   *Note: this feature is only supported currently for the solver SDPT3,
%    and should be considered experimental.*
%
%   CVX_WARMSTART ON     instructs CVX to begin performing warm starts. The
%                        first model solved after this command is issued
%                        will still be "cold" started, but its data will
%                        be saved for the next model.
%   CVX_WARMSTART OFF    instructs CVX to stop performing warm starts.
%   CVX_WARMSTART CLEAR  instructs CVX to clear any current warm start data
%                        but to preserve the current warm start state.
%   CVX_WARMSTART        with no arguments prints the current state.
%   Y = CVX_WARMSTART    returns 'off', 'cold', or 'warm'.
%
%   "Warm starting" is the practice of using a previous model's solution
%   to help initialze the solution process for a new model. With *some*
%   solvers, in *certain* scenarios, this can speed up the solution of
%   the second model. For instance, building a tradeoff curve often 
%   involves solving a sequence of models that are identical except for 
%   a single parameter. Such a scenario is a good candidate for warm start.
%
%   CVX implements warm starting by saving the "raw" solution output of the
%   current solver, and then passing that same data into the solver shim
%   when the next model is ready to be solved. The shim is then responsible
%   for extracting any useful information from this raw solution data.
%
%   It is important to emphasize several points:
%   --- Each solver (or solver shim) is responsible for implementing its
%       own warm start functionality. If a solver does not have warm start
%       implemented, or it chooses not to use the supplied solution data,
%       no warning or error will be issued.
%   --- There is no guarantee that this will actually improve performance.
%       In fact, in some cases, performance is degraded by the use of
%       warm starting. Building an effective warm starting approach for
%       a *general purpose* solver is by no means straightforward.
%   --- If you attempt to use CVX_WARMSTART with two incompatible models,
%       the resulting behavior is undefined. The solver shim may detect
%       this situation and revert to cold start; it may cause an error;
%       or it may cause a solver failure. Only use CVX_WARMSTART with
%       models that are known to be identical in structure.
%   --- Warm start data is saved only for Solved models. Data is not saved
%       from models returning a status of Infeasible, Unbounded, or Failed.
%   --- You may not supply your own initial points; only the immediately
%       previous solution may be used.
%
%   The only solver that currently allows warm starting is SDPT3. If 
%   support solvers add this capability, we will update the corresponding
%   shims to take advantage of it.
%
%   CVX_WARMSTART is a global setting, and takes effect immediately. For
%   instance, it can be called immediately before CVX_END, and it will
%   still apply in full for that model.

global cvx___
cvx_global
if nargin == 0
    if ~cvx___.warmstart{1}, res = 'off'; %#ok
    elseif numel(cvx___.warmstart) == 1, res = 'cold';
    else res = 'warm'; end
    if nargout == 0
        switch res,
            case 'off',  fprintf( 'Warm start is OFF.\n' );
            case 'cold', fprintf( 'Warm start is ON; the next model will be COLD started.\n' );
            case 'warm', fprintf( 'Warm start is ON; the next model will be WARM started.\n' );
        end
        clear res
    end
elseif ~ischar(flag) || size(flag,1) ~= 1,
    error( 'Argument must be a string.' );
else
    switch lower(flag),
        case 'off',   cvx___.warmstart = {false};
        case 'on',    cvx___.warmstart{1} = true;
        case 'clear', cvx___.warmstart(2:end) = [];
        otherwise, error( 'Invalid flag: %s\n', flag );
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
