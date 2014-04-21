function clear( prob, warn )

global cvx___
[ p, pstr ] = verify( prob );
if nargin == 2 && warn && ( ~isempty(pstr.objective) || ~isempty(pstr.variables) || ~isempty(pstr.duals) || nnz(pstr.t_variable) > 1 );
    warning( 'CVX:Empty', 'A non-empty cvx problem already exists in this scope.\n   It is being overwritten.', 1 ); %#ok
    evalin( 'caller', 'cvx_pop( cvx_problem, true )' );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
