function cvx_end

%CVX_END  Completes a cvx specification.
%   CVX_BEGIN marks the end of a new cvx model, and instructs cvx to
%   complete its processing. For standard, complete models, cvx will send
%   a transformed version of the problem to a solver to obtain numeric
%   results, and replace references to cvx variables with numeric values.

try
    evalin( 'caller', 'cvx_verify' );
    evalin( 'caller', 'cvx_finish' );
    evalin( 'caller', 'cvx_pop' );
catch exc
    if ~isequal( exc.identifier, 'CVX:IncompatibleSolver' ),
        evalin( 'caller', 'cvx_pop' );
    end
    rethrow( exc );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
