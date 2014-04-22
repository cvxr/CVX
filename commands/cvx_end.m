function cvx_end

%CVX_END  Completes a cvx specification.
%   CVX_BEGIN marks the end of a new cvx model, and instructs cvx to
%   complete its processing. For standard, complete models, cvx will send
%   a transformed version of the problem to a solver to obtain numeric
%   results, and replace references to cvx variables with numeric values.

if ~isa( evalin( 'caller', 'cvx_problem', '[]' ), 'cvxprob' ), 
    error( 'CVX:NoProblem', 'No model exists in this scope.' ); 
end
try
    evalin( 'caller', 'finish( cvx_problem )' );
catch exc
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc );
    else rethrow( exc ); end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
