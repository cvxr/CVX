function a = gt( x, y )

if isa(y,'cvxin')||~isa(x,'cvxin')||~x.active,
    error( 'CVX error: improper use of the <in> pseudo-operator.' );
end
b = newcnstr( evalin( 'caller', 'cvx_problem', '[]' ), x.value, y, '==' );
if nargout, a = b; end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
