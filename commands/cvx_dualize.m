function cvx_dualize( flag )

%CVX_DUALIZE    Controls CVX dualization.
%   The CVX_DUALIZE command controls CVX's dualization functionality. CVX
%   can solve either the primal or dual version of your problem. CVX uses
%   a heuristic to determine which choice is the most efficient.
%
%   In some instances, however, CVX's automatic determination may not be
%   the best. If this is the case, you may wish to override CVX's default
%   behavior by forcing it to dualize, or forcing it not to dualize. This
%   command gives you that capability:
%      CVX_DUALIZE ON   forces CVX to pass the dual problem to the solver
%      CVX_DUALIZE OFF  forces CVX to pass the primal problem to the solver
%
%   Note that if you make the "wrong" choice, you may find that the model
%   passed to the solver is takes much longer to solve, or fails to achieve
%   the desired accuracy. This is an advanced setting meant to be tried
%   only if you suspect that CVX's automatic determination is problematic
%   for a specific model. It may also be useful as a debugging tool for 
%   solver developers. Most users should never need to use it.
%
%   Some solvers do not accept the dual problem; and models with integer
%   variables do not have a dual. In either of these cases, a call to
%   CVX_DUALIZE ON fails *silently*. That is, CVX ignores the setting and
%   passes the primal problem to the solver without warning.
%
%   Unlike similar commands like CVX_PRECISION, this command can only be
%   supplied within a model---that is, between cvx_begin and cvx_end---and 
%   its behavior applies only to that model.

if nargin < 1 || ~ischar(flag) || size(flag,1) ~= 1 || ~any( strcmpi( flag, { 'on', 'off' } ) )
    error( 'Expecting exactly one argument: ON or OFF.' );
end

global cvx___

try

    evalin( 'caller', 'cvx_verify' );
    cvx___.problems(end).dualize = 1 - 2 * strcmpi( flag, 'off' );
        
catch exc
    
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
    
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
